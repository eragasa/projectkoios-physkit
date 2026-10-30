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
    """Transform constructor-bound physical samples with canonical Numerics."""

    __slots__ = ("_source",)

    def __init__(self, *, source: FinitePeriodicMatrixSamples) -> None:
        if type(source) is not FinitePeriodicMatrixSamples:
            raise TypeError("source must be FinitePeriodicMatrixSamples")
        self._source = source

    def result(self) -> FinitePeriodicMatrixFourierTransformResult:
        """Transform the bound samples without selecting support or tolerance."""
        source = self._source
        numerical = NumpyCenteredUniformMeshMatrixFourierTransformer(
            definition=source._numerical_definition()
        ).result()
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
    """Evaluate one constructor-bound physical matrix series and tensor grid."""

    __slots__ = ("_definition", "_series")

    def __init__(
        self,
        *,
        series: FinitePeriodicMatrixFourierSeries,
        coordinate_axes: tuple[FloatArray, ...],
        maximum_working_elements: int = 1_000_000,
    ) -> None:
        if type(series) is not FinitePeriodicMatrixFourierSeries:
            raise TypeError("series must be FinitePeriodicMatrixFourierSeries")
        self._series = series
        self._definition = (
            NumpyTensorProductFourierMatrixEvaluationDefinition.from_arrays(
                series.source_mesh,
                series._numerical_series(),
                coordinate_axes,
                maximum_working_elements,
            )
        )

    def result(self) -> FinitePeriodicMatrixFourierEvaluation:
        """Evaluate the bound complete or explicitly selected series."""
        series = self._series
        definition = self._definition
        numerical = NumpyTensorProductFourierMatrixSeriesEvaluator(
            definition=definition
        ).result()
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
