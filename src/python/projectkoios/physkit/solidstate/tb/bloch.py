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
    BlockHoppingModel,
    ScalarHoppingModel,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexVectorQuantity,
    MatrixQuantity,
    Unitless,
)


def _validate_fractional_wave_vectors(
    model: ScalarHoppingModel | BlockHoppingModel,
    fractional_wave_vectors: MatrixQuantity,
) -> None:
    """Validate unit and shape at the physics-semantic boundary."""
    if not isinstance(fractional_wave_vectors, MatrixQuantity):
        raise TypeError("fractional_wave_vectors must be MatrixQuantity")
    if not isinstance(fractional_wave_vectors.unit, Unitless):
        raise ValueError("fractional_wave_vectors must be unitless")
    expected_width = model.dimension.value
    if fractional_wave_vectors.magnitude.shape[1] != expected_width:
        raise ValueError(
            f"fractional_wave_vectors must have shape (point_count, {expected_width})"
        )


@dataclass(frozen=True, slots=True, eq=False)
class ScalarHoppingBlochEvaluation(ResultsObject):
    """Retain a hopping model, sampled wave vectors, and correlated Bloch values.

    ``bloch_values`` carries the model's energy unit and may be complex because a
    :class:`ScalarHoppingModel` does not impose Hermiticity. Hermitian-model
    acceptance remains a separate diagnostic or consumer policy.
    """

    model: ScalarHoppingModel
    fractional_wave_vectors: MatrixQuantity
    bloch_values: ComplexVectorQuantity

    def __post_init__(self) -> None:
        """Validate model ownership, coordinate semantics, unit, and shape."""
        if type(self.model) is not ScalarHoppingModel:
            raise TypeError("model must be ScalarHoppingModel")
        _validate_fractional_wave_vectors(self.model, self.fractional_wave_vectors)
        if not isinstance(self.bloch_values, ComplexVectorQuantity):
            raise TypeError("bloch_values must be ComplexVectorQuantity")
        if self.bloch_values.unit != self.model.energy_unit:
            raise ValueError("bloch_values must use the model energy unit")
        point_count = self.fractional_wave_vectors.magnitude.shape[0]
        if self.bloch_values.magnitude.shape != (point_count,):
            raise ValueError("bloch_values must match the wave-vector count")


@dataclass(frozen=True, slots=True, eq=False)
class BlockHoppingBlochEvaluation(ResultsObject):
    """Retain a block model, sampled wave vectors, and Bloch matrices."""

    model: BlockHoppingModel
    fractional_wave_vectors: MatrixQuantity
    bloch_matrices: tuple[ComplexMatrixQuantity, ...]

    def __post_init__(self) -> None:
        """Validate model ownership, units, matrix shapes, and point correlation."""
        if type(self.model) is not BlockHoppingModel:
            raise TypeError("model must be BlockHoppingModel")
        _validate_fractional_wave_vectors(self.model, self.fractional_wave_vectors)
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
    """Adapt a scalar hopping model to the generic finite-series evaluator.

    Parameters
    ----------
    finite_series_evaluator
        Nominal Numerics evaluator to use. When omitted, the reviewed ordered
        NumPy backend is constructed lazily through the Numerics factory.
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
        """Evaluate the scalar Bloch Hamiltonian at supplied fractional vectors.

        Rows of ``fractional_wave_vectors`` are reciprocal coordinates in turns
        relative to the primitive direct-lattice axes. Empty point sets are
        permitted and return an empty energy-unit-bearing complex vector.
        """
        if type(model) is not ScalarHoppingModel:
            raise TypeError("model must be ScalarHoppingModel")
        _validate_fractional_wave_vectors(model, fractional_wave_vectors)

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


class BlockHoppingBlochEvaluator:
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
        model: BlockHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
    ) -> BlockHoppingBlochEvaluation:
        """Evaluate one block model under its Bloch convention.

        Parameters
        ----------
        model : BlockHoppingModel
            Ordered matrix-valued hopping model.
        fractional_wave_vectors : MatrixQuantity
            Unitless reciprocal coordinates with one row per point.

        Returns
        -------
        BlockHoppingBlochEvaluation
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
        if type(model) is not BlockHoppingModel:
            raise TypeError("model must be BlockHoppingModel")
        return self.evaluate_family((model,), fractional_wave_vectors)[0]

    def evaluate_family(
        self,
        models: tuple[BlockHoppingModel, ...],
        fractional_wave_vectors: MatrixQuantity,
    ) -> tuple[BlockHoppingBlochEvaluation, ...]:
        """Evaluate a compatible ordered model family on one coordinate grid.

        Parameters
        ----------
        models : tuple of BlockHoppingModel
            Nonempty ordered models sharing dimension, basis, energy unit, and energy
            reference. Displacement inventories may differ.
        fractional_wave_vectors : MatrixQuantity
            Exact shared unitless coordinate record retained by every result.

        Returns
        -------
        tuple of BlockHoppingBlochEvaluation
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
            type(model) is not BlockHoppingModel for model in models
        ):
            raise TypeError("models must be a tuple of BlockHoppingModel values")
        if not models:
            raise ValueError("models must be nonempty")
        reference = models[0]
        _validate_fractional_wave_vectors(reference, fractional_wave_vectors)
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
        model: BlockHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
    ) -> BlockHoppingBlochEvaluation:
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
        models: tuple[BlockHoppingModel, ...],
        fractional_wave_vectors: MatrixQuantity,
    ) -> tuple[BlockHoppingBlochEvaluation, ...]:
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
        model: BlockHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
        values: np.ndarray,
    ) -> BlockHoppingBlochEvaluation:
        """Restore units and point correlation for one model's matrix stack."""
        return BlockHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=fractional_wave_vectors,
            bloch_matrices=tuple(
                ComplexMatrixQuantity(matrix, model.energy_unit) for matrix in values
            ),
        )
