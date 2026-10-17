"""Morse potential parameter record."""

from dataclasses import dataclass

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
)

_METRE = PhysicalUnit(expression="meter")
_JOULE = PhysicalUnit(expression="joule")
_INVERSE_METRE = PhysicalUnit(expression="1 / meter")


@dataclass(frozen=True, slots=True, kw_only=True)
class MorsePotentialParameters(DataObject):
    """Represent positive Morse energy, inverse-range, and distance parameters."""

    dissociation_energy: ScalarQuantity
    inverse_range: ScalarQuantity
    equilibrium_distance: ScalarQuantity

    def __post_init__(self) -> None:
        """Check every Morse parameter argument."""
        self._check_arg_dissociation_energy()
        self._check_arg_inverse_range()
        self._check_arg_equilibrium_distance()

    def _check_arg_dissociation_energy(self) -> None:
        """Require a positive physical energy."""
        if not isinstance(self.dissociation_energy, ScalarQuantity):
            raise TypeError("dissociation_energy must be ScalarQuantity")
        if not isinstance(self.dissociation_energy.unit, PhysicalUnit):
            raise TypeError("dissociation_energy must use a physical energy unit")
        if not ModelSystemUnitConverter.compatible(
            self.dissociation_energy.unit,
            _JOULE,
        ):
            raise ValueError(
                "dissociation_energy must use units compatible with energy"
            )
        if self.dissociation_energy.magnitude <= 0.0:
            raise ValueError("dissociation_energy must be positive")

    def _check_arg_inverse_range(self) -> None:
        """Require a positive physical inverse length."""
        if not isinstance(self.inverse_range, ScalarQuantity):
            raise TypeError("inverse_range must be ScalarQuantity")
        if not isinstance(self.inverse_range.unit, PhysicalUnit):
            raise TypeError("inverse_range must use a physical inverse-length unit")
        if not ModelSystemUnitConverter.compatible(
            self.inverse_range.unit,
            _INVERSE_METRE,
        ):
            raise ValueError(
                "inverse_range must use units compatible with inverse length"
            )
        if self.inverse_range.magnitude <= 0.0:
            raise ValueError("inverse_range must be positive")

    def _check_arg_equilibrium_distance(self) -> None:
        """Require a positive physical length."""
        if not isinstance(self.equilibrium_distance, ScalarQuantity):
            raise TypeError("equilibrium_distance must be ScalarQuantity")
        if not isinstance(self.equilibrium_distance.unit, PhysicalUnit):
            raise TypeError("equilibrium_distance must use a physical length unit")
        if not ModelSystemUnitConverter.compatible(
            self.equilibrium_distance.unit,
            _METRE,
        ):
            raise ValueError(
                "equilibrium_distance must use units compatible with length"
            )
        if self.equilibrium_distance.magnitude <= 0.0:
            raise ValueError("equilibrium_distance must be positive")
