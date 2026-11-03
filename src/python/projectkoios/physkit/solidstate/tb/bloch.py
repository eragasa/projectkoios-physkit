r"""Bloch evaluation of scalar and block translation-invariant hoppings.

For hopping coefficients :math:`t_{\mathbf R}` or matrices
:math:`T_{\mathbf R}` and fractional reciprocal coordinates :math:`\mathbf q`,
these adapters evaluate

.. math::

   h(\mathbf q) = \sum_{\mathbf R} t_{\mathbf R}
   \exp\!\left(+2\pi i\,\mathbf q\cdot\mathbf R\right)

or the entrywise matrix analogue. The positive exponent and absence of
normalization are explicit conventions. PhysKit owns hopping-model identity,
fractional reciprocal interpretation, basis identity, and energy units.
``projectkoios.numerics`` owns the unit-free ordered scalar and shared-grid
batched finite Fourier-series calculations.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from projectkoios.numerics.finite_periodic.batched_finite_series import (
    BatchedFiniteFourierSeriesEvaluationDefinition,
    BatchedFiniteFourierSeriesEvaluator,
)
from projectkoios.numerics.finite_periodic.finite_series import (
    FiniteFourierSeriesEvaluationDefinition,
    FiniteFourierSeriesEvaluator,
    FiniteFourierSeriesEvaluatorFactory,
    FiniteFourierSeriesMethod,
)
from projectkoios.numerics.finite_periodic.numpy.batched_finite_series import (
    NumpyBatchedFiniteFourierSeriesEvaluator,
)
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.lattice.hopping import (
    BaseBlockHoppingModel,
    DenseBlockHoppingModel,
    ScalarHoppingModel,
    SparseBlockHoppingModel,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    ComplexVectorQuantity,
    MatrixQuantity,
    Unitless,
)


class _BlochFractionalWaveVectorValidator:
    """Validate reciprocal coordinates against one retained hopping model.

    Parameters
    ----------
    model : ScalarHoppingModel or BaseBlockHoppingModel
        Model supplying the required ordinary-lattice dimension.

    Raises
    ------
    TypeError
        If ``model`` is not a supported hopping-model owner.

    Notes
    -----
    This action owns the shared Bloch-coordinate boundary. It performs no coordinate
    conversion and retains no mutable numerical state.
    """

    __slots__ = ("_expected_width",)

    def __init__(self, model: ScalarHoppingModel | BaseBlockHoppingModel) -> None:
        """Retain the exact reciprocal-coordinate width required by ``model``."""
        if not isinstance(model, ScalarHoppingModel | BaseBlockHoppingModel):
            raise TypeError("model must be a scalar or block hopping model")
        self._expected_width = model.dimension.value

    def validate(self, fractional_wave_vectors: MatrixQuantity) -> None:
        """Validate coordinate ownership, unit, and dimensional compatibility.

        Parameters
        ----------
        fractional_wave_vectors : MatrixQuantity
            Candidate reciprocal-coordinate matrix.

        Raises
        ------
        TypeError
            If coordinates are not a matrix quantity.
        ValueError
            If coordinates are dimensional or have the wrong column count.
        """
        if not isinstance(fractional_wave_vectors, MatrixQuantity):
            raise TypeError("fractional_wave_vectors must be MatrixQuantity")
        if not isinstance(fractional_wave_vectors.unit, Unitless):
            raise ValueError("fractional_wave_vectors must be unitless")
        if fractional_wave_vectors.magnitude.shape[1] != self._expected_width:
            raise ValueError(
                "fractional_wave_vectors must have shape "
                f"(point_count, {self._expected_width})"
            )


@dataclass(frozen=True, slots=True, eq=False)
class ScalarHoppingBlochEvaluation(ResultsObject):
    """Correlate a scalar hopping model with sampled complex Bloch values.

    Parameters
    ----------
    model : ScalarHoppingModel
        Exact ordinary-lattice source model.
    fractional_wave_vectors : MatrixQuantity
        Unitless reciprocal coordinates with one row per sampled point.
    bloch_values : ComplexVectorQuantity
        One complex value per coordinate row in the model energy unit.

    Raises
    ------
    TypeError
        If a value lacks its required nominal quantity or model owner.
    ValueError
        If coordinate width, point count, or energy unit is inconsistent with the
        model.

    Notes
    -----
    Values may be complex because the model does not require Hermitian displacement
    pairing. Hermitian-model acceptance remains a separate diagnostic policy.
    """

    model: ScalarHoppingModel
    fractional_wave_vectors: MatrixQuantity
    bloch_values: ComplexVectorQuantity

    def __post_init__(self) -> None:
        """Validate model ownership, coordinate semantics, unit, and shape."""
        if type(self.model) is not ScalarHoppingModel:
            raise TypeError("model must be ScalarHoppingModel")
        _BlochFractionalWaveVectorValidator(self.model).validate(
            self.fractional_wave_vectors
        )
        if not isinstance(self.bloch_values, ComplexVectorQuantity):
            raise TypeError("bloch_values must be ComplexVectorQuantity")
        if self.bloch_values.unit != self.model.energy_unit:
            raise ValueError("bloch_values must use the model energy unit")
        point_count = self.fractional_wave_vectors.magnitude.shape[0]
        if self.bloch_values.magnitude.shape != (point_count,):
            raise ValueError("bloch_values must match the wave-vector count")


@dataclass(frozen=True, slots=True, eq=False)
class DenseBlockHoppingBlochEvaluation(ResultsObject):
    """Correlate a dense hopping model with sampled dense Bloch matrices.

    Parameters
    ----------
    model : DenseBlockHoppingModel
        Exact ordinary-lattice source model.
    fractional_wave_vectors : MatrixQuantity
        Unitless reciprocal coordinates with one row per sampled point.
    bloch_matrices : tuple[ComplexMatrixQuantity, ...]
        One dense square matrix per coordinate row, in matching order.

    Raises
    ------
    TypeError
        If a value lacks its required exact nominal owner.
    ValueError
        If coordinate width, point count, matrix shape, or energy unit is
        inconsistent with the model.

    Notes
    -----
    Structural correlation does not establish Hermiticity, diagonalization, or
    reciprocal-space convergence.
    """

    model: DenseBlockHoppingModel
    fractional_wave_vectors: MatrixQuantity
    bloch_matrices: tuple[ComplexMatrixQuantity, ...]

    def __post_init__(self) -> None:
        """Validate model ownership, units, matrix shapes, and point correlation."""
        if type(self.model) is not DenseBlockHoppingModel:
            raise TypeError("model must be DenseBlockHoppingModel")
        _BlochFractionalWaveVectorValidator(self.model).validate(
            self.fractional_wave_vectors
        )
        if type(self.bloch_matrices) is not tuple or any(
            type(matrix) is not ComplexMatrixQuantity for matrix in self.bloch_matrices
        ):
            raise TypeError(
                "bloch_matrices must be a tuple of ComplexMatrixQuantity values"
            )
        point_count = self.fractional_wave_vectors.magnitude.shape[0]
        if len(self.bloch_matrices) != point_count:
            raise ValueError("bloch_matrices must match the wave-vector count")
        expected_shape = (self.model.basis_size, self.model.basis_size)
        if any(
            matrix.magnitude.shape != expected_shape for matrix in self.bloch_matrices
        ):
            raise ValueError("Bloch matrix shapes must match model basis_size")
        if any(matrix.unit != self.model.energy_unit for matrix in self.bloch_matrices):
            raise ValueError("Bloch matrix units must match the model energy unit")


class ScalarHoppingBlochEvaluator:
    """Adapt a scalar hopping model to a finite Fourier-series evaluator.

    Parameters
    ----------
    finite_series_evaluator : FiniteFourierSeriesEvaluator or None, optional
        Nominal Numerics evaluator. When omitted, the reviewed ordered NumPy backend
        is constructed through the Numerics factory.

    Raises
    ------
    TypeError
        If an injected evaluator lacks the required nominal interface.

    Notes
    -----
    The adapter preserves displacement order, coordinate order, and the model energy
    unit. It does not impose Hermiticity or project values onto the real axis.
    """

    __slots__ = ("_finite_series_evaluator",)

    def __init__(
        self,
        finite_series_evaluator: FiniteFourierSeriesEvaluator | None = None,
    ) -> None:
        """Retain an injected evaluator or construct the built-in NumPy backend."""
        if finite_series_evaluator is None:
            finite_series_evaluator = FiniteFourierSeriesEvaluatorFactory.create(
                FiniteFourierSeriesMethod.NUMPY_DENSE
            )
        elif not isinstance(
            finite_series_evaluator,
            FiniteFourierSeriesEvaluator,
        ):
            raise TypeError(
                "finite_series_evaluator must be FiniteFourierSeriesEvaluator or None"
            )
        self._finite_series_evaluator = finite_series_evaluator

    @property
    def finite_series_evaluator(self) -> FiniteFourierSeriesEvaluator:
        """Return the nominal unit-free evaluator used by this adapter."""
        return self._finite_series_evaluator

    def evaluate(
        self,
        model: ScalarHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
    ) -> ScalarHoppingBlochEvaluation:
        """Evaluate one scalar model at fractional reciprocal coordinates.

        Parameters
        ----------
        model : ScalarHoppingModel
            Ordered scalar hopping model.
        fractional_wave_vectors : MatrixQuantity
            Unitless coordinates of shape ``(point_count, model.dimension.value)``.

        Returns
        -------
        ScalarHoppingBlochEvaluation
            Complex values in coordinate-row order and the model energy unit.

        Raises
        ------
        TypeError
            If the model or coordinates lack their required nominal owners.
        ValueError
            If coordinate units or shape are incompatible with the model.
        FloatingPointError
            If delegated finite-Fourier arithmetic produces nonfinite values.

        Notes
        -----
        Coordinate components are turns relative to primitive reciprocal axes. Empty
        point sets produce an empty energy-unit-bearing complex vector.
        """
        if type(model) is not ScalarHoppingModel:
            raise TypeError("model must be ScalarHoppingModel")
        _BlochFractionalWaveVectorValidator(model).validate(
            fractional_wave_vectors
        )

        # Strip physical wrappers only at the Numerics boundary while preserving
        # canonical displacement and coordinate order exactly.
        definition = FiniteFourierSeriesEvaluationDefinition(
            dimension=model.dimension.value,
            coordinates=tuple(
                tuple(float(component) for component in wave_vector)
                for wave_vector in fractional_wave_vectors.magnitude
            ),
            integer_modes=tuple(term.displacement.components for term in model.terms),
            coefficients=tuple(term.value for term in model.terms),
        )
        finite_series_result = self._finite_series_evaluator.evaluate(definition)

        # Restore the source energy convention and retain the exact coordinate record
        # so result correlation remains explicit and immutable.
        return ScalarHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=fractional_wave_vectors,
            bloch_values=ComplexVectorQuantity(
                magnitude=np.asarray(
                    finite_series_result.values,
                    dtype=np.complex128,
                ),
                unit=model.energy_unit,
            ),
        )


class DenseBlockHoppingBlochEvaluator:
    """Adapt block hopping models to finite Fourier-series owners.

    Parameters
    ----------
    finite_series_evaluator : FiniteFourierSeriesEvaluator or None, optional
        Ordered scalar evaluator retained for compatibility and exact scalar
        reduction. Supplying it without a batched evaluator selects the scalar path.
    batched_finite_series_evaluator : BatchedFiniteFourierSeriesEvaluator or None
        Shared-grid channel evaluator. When both arguments are omitted, the dense
        NumPy batched evaluator is selected for block calculations.

    Notes
    -----
    The batched path shares one point-mode phase matrix across matrix entries and
    compatible model families. Its dense mode reduction can differ from the scalar
    ordered evaluator by binary64 rounding. No diagonalization, Hermiticity decision,
    truncation policy, or basis transformation is performed.
    """

    __slots__ = ("_batched_finite_series_evaluator", "_finite_series_evaluator")

    def __init__(
        self,
        finite_series_evaluator: FiniteFourierSeriesEvaluator | None = None,
        *,
        batched_finite_series_evaluator: (
            BatchedFiniteFourierSeriesEvaluator | None
        ) = None,
    ) -> None:
        """Retain nominal owners and select the default accelerated block path."""
        use_default_batch = finite_series_evaluator is None
        if finite_series_evaluator is None:
            finite_series_evaluator = FiniteFourierSeriesEvaluatorFactory.create(
                FiniteFourierSeriesMethod.NUMPY_DENSE
            )
        elif not isinstance(finite_series_evaluator, FiniteFourierSeriesEvaluator):
            raise TypeError(
                "finite_series_evaluator must be FiniteFourierSeriesEvaluator or None"
            )
        if batched_finite_series_evaluator is None and use_default_batch:
            batched_finite_series_evaluator = NumpyBatchedFiniteFourierSeriesEvaluator()
        elif batched_finite_series_evaluator is not None and not isinstance(
            batched_finite_series_evaluator,
            BatchedFiniteFourierSeriesEvaluator,
        ):
            raise TypeError(
                "batched_finite_series_evaluator must be "
                "BatchedFiniteFourierSeriesEvaluator or None"
            )
        self._finite_series_evaluator = finite_series_evaluator
        self._batched_finite_series_evaluator = batched_finite_series_evaluator

    @property
    def finite_series_evaluator(self) -> FiniteFourierSeriesEvaluator:
        """Return the ordered scalar evaluator retained by this adapter."""
        return self._finite_series_evaluator

    @property
    def batched_finite_series_evaluator(
        self,
    ) -> BatchedFiniteFourierSeriesEvaluator | None:
        """Return the selected shared-grid evaluator, if present."""
        return self._batched_finite_series_evaluator

    def evaluate(
        self,
        model: DenseBlockHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
    ) -> DenseBlockHoppingBlochEvaluation:
        """Evaluate one block model under its Bloch convention.

        Parameters
        ----------
        model : DenseBlockHoppingModel
            Ordered matrix-valued hopping model.
        fractional_wave_vectors : MatrixQuantity
            Unitless reciprocal coordinates with one row per point.

        Returns
        -------
        DenseBlockHoppingBlochEvaluation
            Correlated matrices in coordinate-row order.

        Raises
        ------
        TypeError
            If the model or coordinates have unsupported nominal types.
        ValueError
            If coordinate units or shape are incompatible with the model.
        FloatingPointError
            If delegated finite-Fourier arithmetic becomes nonfinite.

        Notes
        -----
        The default batched backend uses dense channel reduction and need not be
        bitwise identical to the ordered scalar compatibility path.
        """
        if type(model) is not DenseBlockHoppingModel:
            raise TypeError("model must be DenseBlockHoppingModel")
        return self.evaluate_family((model,), fractional_wave_vectors)[0]

    def evaluate_family(
        self,
        models: tuple[DenseBlockHoppingModel, ...],
        fractional_wave_vectors: MatrixQuantity,
    ) -> tuple[DenseBlockHoppingBlochEvaluation, ...]:
        """Evaluate a compatible ordered model family on one coordinate grid.

        Parameters
        ----------
        models : tuple of DenseBlockHoppingModel
            Nonempty ordered models sharing dimension, basis, energy unit, and energy
            reference. Displacement inventories may differ.
        fractional_wave_vectors : MatrixQuantity
            Exact shared unitless coordinate record retained by every result.

        Returns
        -------
        tuple of DenseBlockHoppingBlochEvaluation
            One result per model in supplied order.

        Raises
        ------
        TypeError
            If the family or coordinates have unsupported nominal types.
        ValueError
            If the family is empty, coordinates are invalid, or model metadata is
            incompatible.
        FloatingPointError
            If delegated finite-Fourier arithmetic becomes nonfinite.

        Notes
        -----
        The batched path uses the sorted union of model displacements and exact complex
        zeros for absent coefficients. Empty coordinate grids retain the complete
        ordered result family without invoking numerical backend work.
        """
        if type(models) is not tuple or any(
            type(model) is not DenseBlockHoppingModel for model in models
        ):
            raise TypeError("models must be a tuple of DenseBlockHoppingModel values")
        if not models:
            raise ValueError("models must be nonempty")
        reference = models[0]
        _BlochFractionalWaveVectorValidator(reference).validate(
            fractional_wave_vectors
        )
        if any(model.dimension is not reference.dimension for model in models[1:]):
            raise ValueError("family model dimensions must match")
        if any(model.basis_size != reference.basis_size for model in models[1:]):
            raise ValueError("family model basis sizes must match")
        if any(
            model.basis_identifier != reference.basis_identifier for model in models[1:]
        ):
            raise ValueError("family model basis identifiers must match")
        if any(model.energy_unit != reference.energy_unit for model in models[1:]):
            raise ValueError("family model energy units must match")
        if any(
            model.energy_reference != reference.energy_reference for model in models[1:]
        ):
            raise ValueError("family model energy references must match")
        if self._batched_finite_series_evaluator is None:
            return tuple(
                self._evaluate_ordered(model, fractional_wave_vectors)
                for model in models
            )
        return self._evaluate_batched(models, fractional_wave_vectors)

    def _evaluate_ordered(
        self,
        model: DenseBlockHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
    ) -> DenseBlockHoppingBlochEvaluation:
        """Evaluate one model through the exact ordered scalar owner."""
        coordinates = tuple(
            tuple(float(component) for component in wave_vector)
            for wave_vector in fractional_wave_vectors.magnitude
        )
        integer_modes = tuple(term.displacement.components for term in model.terms)
        values = np.empty(
            (len(coordinates), model.basis_size, model.basis_size),
            dtype=np.complex128,
        )
        for row in range(model.basis_size):
            for column in range(model.basis_size):
                definition = FiniteFourierSeriesEvaluationDefinition(
                    dimension=model.dimension.value,
                    coordinates=coordinates,
                    integer_modes=integer_modes,
                    coefficients=tuple(
                        complex(term.coefficient.magnitude[row, column])
                        for term in model.terms
                    ),
                )
                result = self._finite_series_evaluator.evaluate(definition)
                values[:, row, column] = np.asarray(
                    result.values,
                    dtype=np.complex128,
                )
        return self._result(model, fractional_wave_vectors, values)

    def _evaluate_batched(
        self,
        models: tuple[DenseBlockHoppingModel, ...],
        fractional_wave_vectors: MatrixQuantity,
    ) -> tuple[DenseBlockHoppingBlochEvaluation, ...]:
        """Evaluate compatible models through one shared-grid owner call."""
        reference = models[0]
        coordinates = tuple(
            tuple(float(component) for component in wave_vector)
            for wave_vector in fractional_wave_vectors.magnitude
        )
        # A source-plus-truncations family keeps source displacement order here because
        # model terms and the union are both canonical lexicographic inventories.
        integer_modes = tuple(
            sorted(
                {
                    term.displacement.components
                    for model in models
                    for term in model.terms
                }
            )
        )
        channels: list[tuple[complex, ...]] = []
        for model in models:
            terms = {term.displacement.components: term for term in model.terms}
            for row in range(reference.basis_size):
                for column in range(reference.basis_size):
                    # Exact zeros represent modes omitted from one family member while
                    # retaining one shared phase matrix across the complete family.
                    channels.append(
                        tuple(
                            0.0 + 0.0j
                            if (term := terms.get(mode)) is None
                            else complex(term.coefficient.magnitude[row, column])
                            for mode in integer_modes
                        )
                    )
        definition = BatchedFiniteFourierSeriesEvaluationDefinition(
            dimension=reference.dimension.value,
            coordinates=coordinates,
            integer_modes=integer_modes,
            coefficient_channels=tuple(channels),
        )
        assert self._batched_finite_series_evaluator is not None
        result = self._batched_finite_series_evaluator.evaluate(definition)
        # Owner values are P x (family*n*n). Reshaping preserves coefficient-channel
        # order before the model axis is moved ahead of the reciprocal-point axis.
        values = np.asarray(result.values, dtype=np.complex128).reshape(
            len(coordinates),
            len(models),
            reference.basis_size,
            reference.basis_size,
        )
        return tuple(
            self._result(model, fractional_wave_vectors, values[:, index])
            for index, model in enumerate(models)
        )

    @staticmethod
    def _result(
        model: DenseBlockHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
        values: np.ndarray,
    ) -> DenseBlockHoppingBlochEvaluation:
        """Restore units and point correlation for one model's matrix stack."""
        return DenseBlockHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=fractional_wave_vectors,
            bloch_matrices=tuple(
                ComplexMatrixQuantity(matrix, model.energy_unit) for matrix in values
            ),
        )


@dataclass(frozen=True, slots=True, eq=False)
class SparseBlockHoppingBlochEvaluation(ResultsObject):
    """Correlate a sparse hopping model with sampled sparse Bloch matrices.

    Parameters
    ----------
    model : SparseBlockHoppingModel
        Exact ordinary-lattice source model.
    fractional_wave_vectors : MatrixQuantity
        Unitless reciprocal coordinates with one row per sampled point.
    bloch_matrices : tuple[ComplexSparseMatrixQuantity, ...]
        One canonical sparse square matrix per coordinate row.

    Raises
    ------
    TypeError
        If a value lacks its required exact nominal owner.
    ValueError
        If coordinate width, point count, matrix shape, or energy unit is
        inconsistent with the model.

    Notes
    -----
    The result preserves sparse storage and never materializes a dense matrix.
    Structural correlation does not establish Hermiticity or convergence.
    """

    model: SparseBlockHoppingModel
    fractional_wave_vectors: MatrixQuantity
    bloch_matrices: tuple[ComplexSparseMatrixQuantity, ...]

    def __post_init__(self) -> None:
        """Validate sparse ownership, units, shapes, and point correlation."""
        if type(self.model) is not SparseBlockHoppingModel:
            raise TypeError("model must be SparseBlockHoppingModel")
        _BlochFractionalWaveVectorValidator(self.model).validate(
            self.fractional_wave_vectors
        )
        if type(self.bloch_matrices) is not tuple or any(
            type(matrix) is not ComplexSparseMatrixQuantity
            for matrix in self.bloch_matrices
        ):
            raise TypeError(
                "bloch_matrices must contain ComplexSparseMatrixQuantity values"
            )
        point_count = self.fractional_wave_vectors.magnitude.shape[0]
        if len(self.bloch_matrices) != point_count:
            raise ValueError("bloch_matrices must match the wave-vector count")
        expected_shape = (self.model.basis_size, self.model.basis_size)
        if any(matrix.shape != expected_shape for matrix in self.bloch_matrices):
            raise ValueError("Bloch matrix shapes must match model basis_size")
        if any(matrix.unit != self.model.energy_unit for matrix in self.bloch_matrices):
            raise ValueError("Bloch matrix units must match the model energy unit")


class SparseBlockHoppingBlochEvaluator:
    r"""Evaluate sparse hopping blocks directly in canonical CSR storage.

    Notes
    -----
    For fractional reciprocal coordinate :math:`\mathbf q`, the evaluator computes

    .. math::

       H(\mathbf q) = \sum_{\mathbf R} T(\mathbf R)
       \exp(+2\pi i\,\mathbf q\cdot\mathbf R).

    The positive exponent and absence of normalization match the dense and scalar
    PhysKit Bloch conventions. No dense fallback is permitted.
    """

    __slots__ = ()

    def evaluate(
        self,
        model: SparseBlockHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
    ) -> SparseBlockHoppingBlochEvaluation:
        """Evaluate one sparse model on an ordered reciprocal-coordinate grid.

        Parameters
        ----------
        model : SparseBlockHoppingModel
            Ordered sparse-block hopping model.
        fractional_wave_vectors : MatrixQuantity
            Unitless coordinates of shape ``(point_count, model.dimension.value)``.

        Returns
        -------
        SparseBlockHoppingBlochEvaluation
            Correlated canonical sparse matrices in coordinate-row order.

        Raises
        ------
        TypeError
            If the model or coordinates lack their required nominal owners.
        ValueError
            If coordinate units or shape are incompatible with the model.
        FloatingPointError
            If phase-weighted sparse accumulation produces nonfinite values.
        """
        if type(model) is not SparseBlockHoppingModel:
            raise TypeError("model must be SparseBlockHoppingModel")
        _BlochFractionalWaveVectorValidator(model).validate(
            fractional_wave_vectors
        )

        shape = (model.basis_size, model.basis_size)
        matrices: list[ComplexSparseMatrixQuantity] = []
        for wave_vector in fractional_wave_vectors.magnitude:
            # Start each reciprocal point in sparse storage. The represented pattern
            # grows only through the exact union of phase-weighted source patterns.
            represented = sparse.csr_array(shape, dtype=np.complex128)
            for term in model.terms:
                phase_argument = float(
                    np.dot(
                        wave_vector,
                        np.asarray(term.displacement.components, dtype=np.float64),
                    )
                )
                phase = np.exp(2.0j * np.pi * phase_argument)
                represented = represented + term.coefficient.to_csr() * phase

            # Sparse addition may create duplicate or exactly cancelled entries;
            # restore the canonical quantity contract before retaining the result.
            represented.sum_duplicates()
            represented.eliminate_zeros()
            represented.sort_indices()
            if not np.all(np.isfinite(represented.data.real)) or not np.all(
                np.isfinite(represented.data.imag)
            ):
                raise FloatingPointError(
                    "sparse Bloch evaluation produced nonfinite values"
                )
            matrices.append(
                ComplexSparseMatrixQuantity.from_csr(
                    represented,
                    model.energy_unit,
                )
            )
        return SparseBlockHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=fractional_wave_vectors,
            bloch_matrices=tuple(matrices),
        )
