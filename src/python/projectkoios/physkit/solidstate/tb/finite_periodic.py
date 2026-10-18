"""Unit-aware physical adapters for finite-periodic matrix Fourier operations.

ProjectKoios Numerics owns centered meshes, Fourier arithmetic, tensor-product
contractions, Parseval reductions, and cyclic Hermiticity pairing. This module owns the
physical correlation that those generic kernels deliberately omit: one explicit energy
unit, energy reference, and ordered-basis identity across samples, coefficients,
evaluations, and diagnostics.

The records and actions choose no support, cutoff, tolerance, or scientific acceptance
policy. Reduced coordinates and integer modes are dimensionless; matrix samples,
coefficients, evaluated matrices, and Hermiticity defects carry the declared energy
unit. Parseval quantities carry its square.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Self, final

import numpy as np
from numpy.typing import NDArray

from projectkoios.numerics.finite_periodic.meshes import CenteredUniformMesh
from projectkoios.numerics.finite_periodic.numpy.matrix_evaluation import (
    NumpyFiniteFourierMatrixEvaluationResult,
)
from projectkoios.numerics.finite_periodic.numpy.matrix_series import (
    NumpyFiniteFourierMatrixSeries,
)
from projectkoios.numerics.finite_periodic.numpy.tensor_product_evaluation import (
    NumpyTensorProductFourierMatrixEvaluationDefinition,
    NumpyTensorProductFourierMatrixSeriesEvaluator,
)
from projectkoios.numerics.finite_periodic.numpy.transform_diagnostics import (
    NumpyCenteredUniformMeshHermiticityAnalyzer,
    NumpyFinitePeriodicParsevalAnalyzer,
)
from projectkoios.numerics.finite_periodic.numpy.transforms import (
    NumpyCenteredUniformMeshMatrixFourierTransformDefinition,
    NumpyCenteredUniformMeshMatrixFourierTransformer,
    NumpyCenteredUniformMeshMatrixFourierTransformResult,
)
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.lattice.hopping import SolidStateUnit
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity, Unitless

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def _identity(value: str, name: str) -> None:
    """Validate one nonempty exact string identity."""
    if type(value) is not str:
        raise TypeError(f"{name} must be a string")
    if not value:
        raise ValueError(f"{name} must be nonempty")


def _unit(value: SolidStateUnit, name: str) -> None:
    """Validate one exact PhysKit solid-state unit identity."""
    if not isinstance(value, PhysicalUnit | Unitless):
        raise TypeError(f"{name} must be PhysicalUnit or Unitless")


def _squared_unit(unit: SolidStateUnit) -> SolidStateUnit:
    """Return the explicit square of one solid-state unit."""
    if isinstance(unit, Unitless):
        return Unitless()
    return PhysicalUnit(f"({unit.expression}) ** 2")


@final
@dataclass(frozen=True, slots=True, eq=False)
class FinitePeriodicMatrixSamples:
    """Represent complete matrix samples with correlated physical identities.

    Parameters
    ----------
    mesh
        Dimensionless centered reduced-coordinate mesh.
    matrix_samples
        Canonical immutable complex128 tensor with shape
        ``(mesh.size,) * mesh.dimension + (n, n)``.
    energy_unit
        Explicit physical or dimensionless unit shared by every matrix entry.
    energy_reference
        Nonempty identity for the common energy zero.
    basis_identifier
        Nonempty identity for the ordered row and column basis.
    """

    mesh: CenteredUniformMesh
    matrix_samples: ComplexArray = field(repr=False)
    energy_unit: SolidStateUnit
    energy_reference: str
    basis_identifier: str

    @classmethod
    def from_array(
        cls,
        *,
        mesh: CenteredUniformMesh,
        matrix_samples: ComplexArray,
        energy_unit: SolidStateUnit,
        energy_reference: str,
        basis_identifier: str,
    ) -> Self:
        """Construct from an exact complex128 array using a defensive copy."""
        numerical = NumpyCenteredUniformMeshMatrixFourierTransformDefinition.from_array(
            mesh,
            matrix_samples,
        )
        return cls(
            mesh=numerical.mesh,
            matrix_samples=numerical.matrix_samples,
            energy_unit=energy_unit,
            energy_reference=energy_reference,
            basis_identifier=basis_identifier,
        )

    def __post_init__(self) -> None:
        """Validate mesh shape, immutable storage, units, and physical identities."""
        if type(self.mesh) is not CenteredUniformMesh:
            raise TypeError("mesh must be CenteredUniformMesh")
        NumpyCenteredUniformMeshMatrixFourierTransformDefinition.from_immutable_array(
            self.mesh,
            self.matrix_samples,
        )
        _unit(self.energy_unit, "energy_unit")
        _identity(self.energy_reference, "energy_reference")
        _identity(self.basis_identifier, "basis_identifier")

    @property
    def matrix_size(self) -> int:
        """Return the common square-matrix size."""
        return int(self.matrix_samples.shape[-1])

    def _numerical_definition(
        self,
    ) -> NumpyCenteredUniformMeshMatrixFourierTransformDefinition:
        """Return the exact unit-free owner definition without copying storage."""
        definition_type = NumpyCenteredUniformMeshMatrixFourierTransformDefinition
        return definition_type.from_immutable_array(
            self.mesh,
            self.matrix_samples,
        )


@final
@dataclass(frozen=True, slots=True, eq=False)
class FinitePeriodicMatrixFourierSeries:
    """Represent selected matrix Fourier coefficients with physical correlation.

    The series may contain a complete centered inventory or an explicitly selected
    ordered subset. Support selection is supplied by a consumer; this record neither
    derives nor recommends a support.
    """

    source_mesh: CenteredUniformMesh
    integer_modes: IntArray = field(repr=False)
    matrix_coefficients: ComplexArray = field(repr=False)
    energy_unit: SolidStateUnit
    energy_reference: str
    basis_identifier: str

    @classmethod
    def from_arrays(
        cls,
        *,
        source_mesh: CenteredUniformMesh,
        integer_modes: IntArray,
        matrix_coefficients: ComplexArray,
        energy_unit: SolidStateUnit,
        energy_reference: str,
        basis_identifier: str,
    ) -> Self:
        """Construct from exact arrays using immutable defensive copies."""
        if type(source_mesh) is not CenteredUniformMesh:
            raise TypeError("source_mesh must be CenteredUniformMesh")
        numerical = NumpyFiniteFourierMatrixSeries.from_arrays(
            source_mesh.dimension,
            integer_modes,
            matrix_coefficients,
        )
        return cls(
            source_mesh=source_mesh,
            integer_modes=numerical.integer_modes,
            matrix_coefficients=numerical.matrix_coefficients,
            energy_unit=energy_unit,
            energy_reference=energy_reference,
            basis_identifier=basis_identifier,
        )

    def __post_init__(self) -> None:
        """Validate series storage, mesh bounds, units, and physical identities."""
        if type(self.source_mesh) is not CenteredUniformMesh:
            raise TypeError("source_mesh must be CenteredUniformMesh")
        numerical = self._numerical_series()
        lower = -(self.source_mesh.size // 2)
        upper = (self.source_mesh.size - 1) // 2
        if bool(
            (
                (numerical.integer_modes < lower) | (numerical.integer_modes > upper)
            ).any()
        ):
            raise ValueError("integer modes must lie inside the source mesh inventory")
        mode_rows = tuple(
            tuple(int(component) for component in row)
            for row in numerical.integer_modes
        )
        if mode_rows != tuple(sorted(set(mode_rows))):
            raise ValueError(
                "integer modes must be lexicographically ordered and unique"
            )
        _unit(self.energy_unit, "energy_unit")
        _identity(self.energy_reference, "energy_reference")
        _identity(self.basis_identifier, "basis_identifier")

    @property
    def dimension(self) -> int:
        """Return the reduced-coordinate dimension."""
        return self.source_mesh.dimension

    @property
    def term_count(self) -> int:
        """Return the ordered coefficient count."""
        return int(self.integer_modes.shape[0])

    @property
    def matrix_size(self) -> int:
        """Return the common square-matrix size."""
        return int(self.matrix_coefficients.shape[-1])

    def take(self, term_indices: tuple[int, ...]) -> Self:
        """Return an explicitly selected ordered subset without choosing policy.

        ``term_indices`` must be strictly increasing and unique, so the selected series
        preserves the source inventory order. The empty selection is supported.
        """
        if type(term_indices) is not tuple or any(
            type(index) is not int for index in term_indices
        ):
            raise TypeError("term_indices must be a tuple of built-in integers")
        if term_indices != tuple(sorted(set(term_indices))):
            raise ValueError("term_indices must be strictly increasing and unique")
        if any(index < 0 or index >= self.term_count for index in term_indices):
            raise ValueError("term_indices must identify represented terms")
        indices = np.asarray(term_indices, dtype=np.int64)
        return FinitePeriodicMatrixFourierSeries.from_arrays(
            source_mesh=self.source_mesh,
            integer_modes=self.integer_modes[indices],
            matrix_coefficients=self.matrix_coefficients[indices],
            energy_unit=self.energy_unit,
            energy_reference=self.energy_reference,
            basis_identifier=self.basis_identifier,
        )

    def _numerical_series(self) -> NumpyFiniteFourierMatrixSeries:
        """Return the exact unit-free owner series without copying storage."""
        return NumpyFiniteFourierMatrixSeries.from_immutable_arrays(
            self.source_mesh.dimension,
            self.integer_modes,
            self.matrix_coefficients,
        )


@final
@dataclass(frozen=True, slots=True, eq=False)
class FinitePeriodicMatrixFourierTransformResult(ResultsObject):
    """Retain a complete transform and its source physical field."""

    source: FinitePeriodicMatrixSamples
    series: FinitePeriodicMatrixFourierSeries
    backend_identifier: str
    backend_version: str

    def __post_init__(self) -> None:
        """Validate complete inventory, physical correlation, and backend evidence."""
        if type(self.source) is not FinitePeriodicMatrixSamples:
            raise TypeError("source must be FinitePeriodicMatrixSamples")
        if type(self.series) is not FinitePeriodicMatrixFourierSeries:
            raise TypeError("series must be FinitePeriodicMatrixFourierSeries")
        if self.series.source_mesh is not self.source.mesh:
            raise ValueError("series source_mesh must be the exact source mesh")
        if self.series.matrix_size != self.source.matrix_size:
            raise ValueError("series and source matrix sizes must agree")
        if (
            self.series.energy_unit != self.source.energy_unit
            or self.series.energy_reference != self.source.energy_reference
            or self.series.basis_identifier != self.source.basis_identifier
        ):
            raise ValueError("series and source physical identities must agree")
        self._numerical_result()

    def _numerical_result(
        self,
    ) -> NumpyCenteredUniformMeshMatrixFourierTransformResult:
        """Return the exact unit-free owner result over maintained storage."""
        return NumpyCenteredUniformMeshMatrixFourierTransformResult(
            mesh=self.source.mesh,
            series=self.series._numerical_series(),
            backend_identifier=self.backend_identifier,
            backend_version=self.backend_version,
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
        return FinitePeriodicMatrixFourierTransformResult(
            source=source,
            series=series,
            backend_identifier=numerical.backend_identifier,
            backend_version=numerical.backend_version,
        )


@final
@dataclass(frozen=True, slots=True, eq=False)
class FinitePeriodicMatrixFourierEvaluation(ResultsObject):
    """Retain tensor-grid matrix values correlated with one physical series."""

    series: FinitePeriodicMatrixFourierSeries
    coordinate_axes: tuple[FloatArray, ...] = field(repr=False)
    matrix_values: ComplexArray = field(repr=False)
    maximum_working_elements: int
    backend_identifier: str
    backend_version: str

    def __post_init__(self) -> None:
        """Validate axes, values, physical correlation, and backend evidence."""
        if type(self.series) is not FinitePeriodicMatrixFourierSeries:
            raise TypeError("series must be FinitePeriodicMatrixFourierSeries")
        definition = (
            NumpyTensorProductFourierMatrixEvaluationDefinition.from_immutable_arrays(
                self.series.source_mesh,
                self.series._numerical_series(),
                self.coordinate_axes,
                self.maximum_working_elements,
            )
        )
        NumpyFiniteFourierMatrixEvaluationResult(
            dimension=self.series.dimension,
            point_shape=definition.point_shape,
            term_count=self.series.term_count,
            matrix_values=self.matrix_values,
            backend_identifier=self.backend_identifier,
            backend_version=self.backend_version,
        )
        expected_matrix_shape = (self.series.matrix_size, self.series.matrix_size)
        if self.matrix_values.shape[-2:] != expected_matrix_shape:
            raise ValueError("matrix value shapes must match the series matrix size")

    @property
    def energy_unit(self) -> SolidStateUnit:
        """Return the exact series energy unit carried by every matrix value."""
        return self.series.energy_unit

    @property
    def energy_reference(self) -> str:
        """Return the exact series energy-reference identity."""
        return self.series.energy_reference

    @property
    def basis_identifier(self) -> str:
        """Return the exact ordered-basis identity."""
        return self.series.basis_identifier

    @property
    def point_shape(self) -> tuple[int, ...]:
        """Return the tensor-product coordinate shape."""
        return tuple(int(axis.size) for axis in self.coordinate_axes)


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
        return FinitePeriodicMatrixFourierEvaluation(
            series=series,
            coordinate_axes=definition.coordinate_axes,
            matrix_values=numerical.matrix_values,
            maximum_working_elements=definition.maximum_working_elements,
            backend_identifier=numerical.backend_identifier,
            backend_version=numerical.backend_version,
        )


@final
@dataclass(frozen=True, slots=True, eq=False)
class FinitePeriodicMatrixParsevalDiagnostics(ResultsObject):
    """Retain unit-aware Parseval evidence for one complete transform."""

    transform: FinitePeriodicMatrixFourierTransformResult
    sample_mean_squared_frobenius_norm: ScalarQuantity
    coefficient_sum_squared_frobenius_norm: ScalarQuantity
    absolute_defect: ScalarQuantity

    def __post_init__(self) -> None:
        """Validate exact canonical values and squared-energy units."""
        if type(self.transform) is not FinitePeriodicMatrixFourierTransformResult:
            raise TypeError(
                "transform must be FinitePeriodicMatrixFourierTransformResult"
            )
        quantities = (
            self.sample_mean_squared_frobenius_norm,
            self.coefficient_sum_squared_frobenius_norm,
            self.absolute_defect,
        )
        if any(type(quantity) is not ScalarQuantity for quantity in quantities):
            raise TypeError("Parseval values must be ScalarQuantity values")
        squared = _squared_unit(self.transform.source.energy_unit)
        if any(quantity.unit != squared for quantity in quantities):
            raise ValueError("Parseval values must use the squared energy unit")
        expected = NumpyFinitePeriodicParsevalAnalyzer().analyze(
            self.transform.source._numerical_definition(),
            self.transform._numerical_result(),
        )
        expected_values = (
            expected.sample_mean_squared_frobenius_norm,
            expected.coefficient_sum_squared_frobenius_norm,
            expected.absolute_defect,
        )
        if tuple(quantity.magnitude for quantity in quantities) != expected_values:
            raise ValueError("Parseval values must equal the canonical diagnostics")


@final
class FinitePeriodicMatrixParsevalAnalyzer:
    """Attach squared-energy units to canonical finite Parseval diagnostics."""

    __slots__ = ()

    def analyze(
        self,
        transform: FinitePeriodicMatrixFourierTransformResult,
    ) -> FinitePeriodicMatrixParsevalDiagnostics:
        """Return threshold-free unit-aware Parseval evidence."""
        if type(transform) is not FinitePeriodicMatrixFourierTransformResult:
            raise TypeError(
                "transform must be FinitePeriodicMatrixFourierTransformResult"
            )
        numerical = NumpyFinitePeriodicParsevalAnalyzer().analyze(
            transform.source._numerical_definition(),
            transform._numerical_result(),
        )
        unit = _squared_unit(transform.source.energy_unit)
        return FinitePeriodicMatrixParsevalDiagnostics(
            transform=transform,
            sample_mean_squared_frobenius_norm=ScalarQuantity(
                numerical.sample_mean_squared_frobenius_norm,
                unit,
            ),
            coefficient_sum_squared_frobenius_norm=ScalarQuantity(
                numerical.coefficient_sum_squared_frobenius_norm,
                unit,
            ),
            absolute_defect=ScalarQuantity(numerical.absolute_defect, unit),
        )


@final
@dataclass(frozen=True, slots=True, eq=False)
class FinitePeriodicMatrixHermiticityDiagnostics(ResultsObject):
    """Retain unit-aware cyclic coefficient-Hermiticity evidence."""

    transform: FinitePeriodicMatrixFourierTransformResult
    twisted_cyclic_maximum_frobenius_defect: ScalarQuantity
    ordinary_maximum_frobenius_defect: ScalarQuantity

    def __post_init__(self) -> None:
        """Validate exact canonical defects and energy units."""
        if type(self.transform) is not FinitePeriodicMatrixFourierTransformResult:
            raise TypeError(
                "transform must be FinitePeriodicMatrixFourierTransformResult"
            )
        quantities = (
            self.twisted_cyclic_maximum_frobenius_defect,
            self.ordinary_maximum_frobenius_defect,
        )
        if any(type(quantity) is not ScalarQuantity for quantity in quantities):
            raise TypeError("Hermiticity defects must be ScalarQuantity values")
        if any(
            quantity.unit != self.transform.source.energy_unit
            for quantity in quantities
        ):
            raise ValueError("Hermiticity defects must use the energy unit")
        expected = NumpyCenteredUniformMeshHermiticityAnalyzer().analyze(
            self.transform._numerical_result()
        )
        expected_values = (
            expected.twisted_cyclic_maximum_frobenius_defect,
            expected.ordinary_maximum_frobenius_defect,
        )
        if tuple(quantity.magnitude for quantity in quantities) != expected_values:
            raise ValueError("Hermiticity defects must equal canonical diagnostics")

    @property
    def has_nyquist_boundary(self) -> bool:
        """Return whether each centered axis has one Nyquist representative."""
        return self.transform.source.mesh.size % 2 == 0

    @property
    def ordinary_comparison_covers_complete_inventory(self) -> bool:
        """Return whether ordinary negation covers every represented mode.

        Odd centered meshes have no Nyquist singleton, so ordinary pairing covers the
        complete inventory. Even meshes require the separately reported twisted cyclic
        comparison at their lower Nyquist boundaries.
        """
        return not self.has_nyquist_boundary


@final
class FinitePeriodicMatrixHermiticityAnalyzer:
    """Attach energy units to canonical cyclic Hermiticity diagnostics."""

    __slots__ = ()

    def analyze(
        self,
        transform: FinitePeriodicMatrixFourierTransformResult,
    ) -> FinitePeriodicMatrixHermiticityDiagnostics:
        """Return threshold-free unit-aware coefficient defects."""
        if type(transform) is not FinitePeriodicMatrixFourierTransformResult:
            raise TypeError(
                "transform must be FinitePeriodicMatrixFourierTransformResult"
            )
        numerical = NumpyCenteredUniformMeshHermiticityAnalyzer().analyze(
            transform._numerical_result()
        )
        unit = transform.source.energy_unit
        return FinitePeriodicMatrixHermiticityDiagnostics(
            transform=transform,
            twisted_cyclic_maximum_frobenius_defect=ScalarQuantity(
                numerical.twisted_cyclic_maximum_frobenius_defect,
                unit,
            ),
            ordinary_maximum_frobenius_defect=ScalarQuantity(
                numerical.ordinary_maximum_frobenius_defect,
                unit,
            ),
        )


__all__ = [
    "FinitePeriodicMatrixFourierEvaluation",
    "FinitePeriodicMatrixFourierSeries",
    "FinitePeriodicMatrixFourierTensorEvaluator",
    "FinitePeriodicMatrixFourierTransformResult",
    "FinitePeriodicMatrixFourierTransformer",
    "FinitePeriodicMatrixHermiticityAnalyzer",
    "FinitePeriodicMatrixHermiticityDiagnostics",
    "FinitePeriodicMatrixParsevalAnalyzer",
    "FinitePeriodicMatrixParsevalDiagnostics",
    "FinitePeriodicMatrixSamples",
]
