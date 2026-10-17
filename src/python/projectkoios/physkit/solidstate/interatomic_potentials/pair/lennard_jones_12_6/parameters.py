"""Lennard--Jones 12-6 parameter record."""

from dataclasses import dataclass

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
)

_METRE = PhysicalUnit(expression="meter")
_JOULE = PhysicalUnit(expression="joule")


@dataclass(frozen=True, slots=True, kw_only=True)
class LennardJonesPotentialParameters(DataObject):
    """Represent positive Lennard--Jones well-depth and distance parameters."""

    well_depth: ScalarQuantity
    zero_crossing_distance: ScalarQuantity

    def __post_init__(self) -> None:
        """Check every Lennard--Jones parameter argument."""
        self._check_arg_well_depth()
        self._check_arg_zero_crossing_distance()

    def _check_arg_well_depth(self) -> None:
        """Require a positive physical energy."""
        if not isinstance(self.well_depth, ScalarQuantity):
            raise TypeError("well_depth must be ScalarQuantity")
        if not isinstance(self.well_depth.unit, PhysicalUnit):
            raise TypeError("well_depth must use a physical energy unit")
        if not ModelSystemUnitConverter.compatible(
            self.well_depth.unit,
            _JOULE,
        ):
            raise ValueError("well_depth must use units compatible with energy")
        if self.well_depth.magnitude <= 0.0:
            raise ValueError("well_depth must be positive")

    def _check_arg_zero_crossing_distance(self) -> None:
        """Require a positive physical length."""
        if not isinstance(self.zero_crossing_distance, ScalarQuantity):
            raise TypeError("zero_crossing_distance must be ScalarQuantity")
        if not isinstance(self.zero_crossing_distance.unit, PhysicalUnit):
            raise TypeError("zero_crossing_distance must use a physical length unit")
        if not ModelSystemUnitConverter.compatible(
            self.zero_crossing_distance.unit,
            _METRE,
        ):
            raise ValueError(
                "zero_crossing_distance must use units compatible with length"
            )
        if self.zero_crossing_distance.magnitude <= 0.0:
            raise ValueError("zero_crossing_distance must be positive")

    @property
    def equilibrium_distance(self) -> ScalarQuantity:
        r"""Return $r_0=2^{1/6}\sigma$ in the input length unit."""
        return ScalarQuantity(
            magnitude=float(2.0 ** (1.0 / 6.0) * self.zero_crossing_distance.magnitude),
            unit=self.zero_crossing_distance.unit,
        )
