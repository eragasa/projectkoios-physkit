"""Unit-aware diagnostics for complete finite-periodic matrix transforms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Self, final

from projectkoios.numerics.finite_periodic.numpy.transform_diagnostics import (
    NumpyCenteredUniformMeshHermiticityAnalyzer,
    NumpyFinitePeriodicParsevalAnalyzer,
)
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ModelSystemUnit,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
)

from .records import FinitePeriodicMatrixFourierTransformResult


def _squared_unit(unit: ModelSystemUnit) -> ModelSystemUnit:
    """Return the explicit square of one solid-state unit."""
    if isinstance(unit, Unitless):
        return Unitless()
    if not isinstance(unit, PhysicalUnit):
        raise TypeError("unit must be PhysicalUnit or Unitless")
    return PhysicalUnit(f"({unit.expression}) ** 2")


@final
@dataclass(frozen=True, slots=True, eq=False, init=False)
class FinitePeriodicMatrixParsevalDiagnostics(ResultsObject):
    """Retain unit-aware Parseval evidence for one complete transform.

    Computed diagnostics are constructed only by
    :class:`FinitePeriodicMatrixParsevalAnalyzer`.
    """

    transform: FinitePeriodicMatrixFourierTransformResult
    sample_mean_squared_frobenius_norm: ScalarQuantity
    coefficient_sum_squared_frobenius_norm: ScalarQuantity
    absolute_defect: ScalarQuantity

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        """Reject direct construction outside the owning diagnostic action."""
        raise TypeError(
            "FinitePeriodicMatrixParsevalDiagnostics must be created by its "
            "owning action"
        )

    @classmethod
    def _from_action(
        cls,
        *,
        transform: FinitePeriodicMatrixFourierTransformResult,
        sample_mean_squared_frobenius_norm: ScalarQuantity,
        coefficient_sum_squared_frobenius_norm: ScalarQuantity,
        absolute_defect: ScalarQuantity,
    ) -> Self:
        """Construct one result for the owning package action."""
        result = object.__new__(cls)
        object.__setattr__(result, "transform", transform)
        object.__setattr__(
            result,
            "sample_mean_squared_frobenius_norm",
            sample_mean_squared_frobenius_norm,
        )
        object.__setattr__(
            result,
            "coefficient_sum_squared_frobenius_norm",
            coefficient_sum_squared_frobenius_norm,
        )
        object.__setattr__(result, "absolute_defect", absolute_defect)
        result._validate()
        return result

    def _validate(self) -> None:
        """Validate the represented relation and squared-energy units."""
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
        expected_defect = abs(
            self.sample_mean_squared_frobenius_norm.magnitude
            - self.coefficient_sum_squared_frobenius_norm.magnitude
        )
        if self.absolute_defect.magnitude != expected_defect:
            raise ValueError("absolute_defect must equal the represented energy defect")


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
        return FinitePeriodicMatrixParsevalDiagnostics._from_action(
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
@dataclass(frozen=True, slots=True, eq=False, init=False)
class FinitePeriodicMatrixHermiticityDiagnostics(ResultsObject):
    """Retain unit-aware cyclic coefficient-Hermiticity evidence.

    Computed diagnostics are constructed only by
    :class:`FinitePeriodicMatrixHermiticityAnalyzer`.
    """

    transform: FinitePeriodicMatrixFourierTransformResult
    twisted_cyclic_maximum_frobenius_defect: ScalarQuantity
    ordinary_maximum_frobenius_defect: ScalarQuantity

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        """Reject direct construction outside the owning diagnostic action."""
        raise TypeError(
            "FinitePeriodicMatrixHermiticityDiagnostics must be created by its "
            "owning action"
        )

    @classmethod
    def _from_action(
        cls,
        *,
        transform: FinitePeriodicMatrixFourierTransformResult,
        twisted_cyclic_maximum_frobenius_defect: ScalarQuantity,
        ordinary_maximum_frobenius_defect: ScalarQuantity,
    ) -> Self:
        """Construct one result for the owning package action."""
        result = object.__new__(cls)
        object.__setattr__(result, "transform", transform)
        object.__setattr__(
            result,
            "twisted_cyclic_maximum_frobenius_defect",
            twisted_cyclic_maximum_frobenius_defect,
        )
        object.__setattr__(
            result,
            "ordinary_maximum_frobenius_defect",
            ordinary_maximum_frobenius_defect,
        )
        result._validate()
        return result

    def _validate(self) -> None:
        """Validate exact energy units and canonical quantity types."""
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
        return FinitePeriodicMatrixHermiticityDiagnostics._from_action(
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
    "FinitePeriodicMatrixHermiticityAnalyzer",
    "FinitePeriodicMatrixHermiticityDiagnostics",
    "FinitePeriodicMatrixParsevalAnalyzer",
    "FinitePeriodicMatrixParsevalDiagnostics",
]
