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
``projectkoios.numerics`` owns the unit-free ordered finite Fourier-series
calculation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.numerics.finite_periodic.finite_series import (
    FiniteFourierSeriesEvaluationDefinition,
    FiniteFourierSeriesEvaluator,
    FiniteFourierSeriesEvaluatorFactory,
    FiniteFourierSeriesMethod,
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
    """Adapt a block hopping model to the scalar finite-series owner.

    Each matrix entry is one ordered finite Fourier series over the shared model
    displacements. No diagonalization, Hermiticity decision, truncation, or basis
    transformation is performed.
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
        elif not isinstance(finite_series_evaluator, FiniteFourierSeriesEvaluator):
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
        model: BlockHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
    ) -> BlockHoppingBlochEvaluation:
        """Evaluate every matrix entry under the model's Bloch convention.

        Rows of ``fractional_wave_vectors`` are unitless reciprocal coordinates in
        turns relative to the primitive direct-lattice axes. An empty point set
        returns an empty tuple of Bloch matrices.
        """
        if type(model) is not BlockHoppingModel:
            raise TypeError("model must be BlockHoppingModel")
        _validate_fractional_wave_vectors(model, fractional_wave_vectors)
        coordinates = tuple(
            tuple(float(component) for component in wave_vector)
            for wave_vector in fractional_wave_vectors.magnitude
        )
        integer_modes = tuple(term.displacement.components for term in model.terms)
        values = np.empty(
            (
                len(coordinates),
                model.basis_size,
                model.basis_size,
            ),
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
        return BlockHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=fractional_wave_vectors,
            bloch_matrices=tuple(
                ComplexMatrixQuantity(matrix, model.energy_unit) for matrix in values
            ),
        )
