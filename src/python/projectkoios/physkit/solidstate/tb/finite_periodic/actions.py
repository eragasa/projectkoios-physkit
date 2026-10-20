"""Unit-aware actions for finite-periodic matrix Fourier operations."""

from __future__ import annotations

from typing import final

from projectkoios.numerics.finite_periodic.numpy.tensor_product_evaluation import (
    NumpyTensorProductFourierMatrixEvaluationDefinition,
    NumpyTensorProductFourierMatrixSeriesEvaluator,
)
from projectkoios.numerics.finite_periodic.numpy.transforms import (
    NumpyCenteredUniformMeshMatrixFourierTransformer,
)

from .records import (
    FinitePeriodicMatrixFourierEvaluation,
    FinitePeriodicMatrixFourierSeries,
    FinitePeriodicMatrixFourierTransformResult,
    FinitePeriodicMatrixSamples,
    FloatArray,
)


@final
class FinitePeriodicMatrixFourierTransformer:
    """Apply the canonical Numerics transform while preserving physical units."""

    __slots__ = ()

    def transform(
        self,
        source: FinitePeriodicMatrixSamples,
    ) -> FinitePeriodicMatrixFourierTransformResult:
        """Transform complete samples without selecting support or tolerance."""
        if type(source) is not FinitePeriodicMatrixSamples:
            raise TypeError("source must be FinitePeriodicMatrixSamples")
        numerical = NumpyCenteredUniformMeshMatrixFourierTransformer().transform(
            source._numerical_definition()
        )
        series = FinitePeriodicMatrixFourierSeries(
            source_mesh=source.mesh,
            integer_modes=numerical.series.integer_modes,
            matrix_coefficients=numerical.series.matrix_coefficients,
            energy_unit=source.energy_unit,
            energy_reference=source.energy_reference,
            basis_identifier=source.basis_identifier,
        )
        return FinitePeriodicMatrixFourierTransformResult._from_action(
            source=source,
            series=series,
            backend_identifier=numerical.backend_identifier,
            backend_version=numerical.backend_version,
        )


@final
class FinitePeriodicMatrixFourierTensorEvaluator:
    """Evaluate a physical matrix series through the canonical tensor kernel."""

    __slots__ = ()

    def evaluate(
        self,
        series: FinitePeriodicMatrixFourierSeries,
        coordinate_axes: tuple[FloatArray, ...],
        maximum_working_elements: int = 1_000_000,
    ) -> FinitePeriodicMatrixFourierEvaluation:
        """Evaluate one complete or explicitly selected series on a tensor grid."""
        if type(series) is not FinitePeriodicMatrixFourierSeries:
            raise TypeError("series must be FinitePeriodicMatrixFourierSeries")
        definition = NumpyTensorProductFourierMatrixEvaluationDefinition.from_arrays(
            series.source_mesh,
            series._numerical_series(),
            coordinate_axes,
            maximum_working_elements,
        )
        numerical = NumpyTensorProductFourierMatrixSeriesEvaluator().evaluate(
            definition
        )
        return FinitePeriodicMatrixFourierEvaluation._from_action(
            series=series,
            coordinate_axes=definition.coordinate_axes,
            matrix_values=numerical.matrix_values,
            maximum_working_elements=definition.maximum_working_elements,
            backend_identifier=numerical.backend_identifier,
            backend_version=numerical.backend_version,
        )


__all__ = [
    "FinitePeriodicMatrixFourierTensorEvaluator",
    "FinitePeriodicMatrixFourierTransformer",
]
