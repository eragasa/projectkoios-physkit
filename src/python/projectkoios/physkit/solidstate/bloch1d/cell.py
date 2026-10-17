r"""Physical one-dimensional primitive cell for Bloch representations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    UnitSystem,
)

_METRE = PhysicalUnit("meter")


@dataclass(frozen=True, slots=True, kw_only=True)
class BlochPrimitiveCell1D(DataObject):
    """Represent a nonzero oriented one-dimensional primitive-cell length."""

    oriented_length: ScalarQuantity
    unit_system: UnitSystem

    def __post_init__(self) -> None:
        """Check every primitive-cell argument."""
        self._check_arg_oriented_length()
        self._check_arg_unit_system()

    def _check_arg_oriented_length(self) -> None:
        if not isinstance(self.oriented_length, ScalarQuantity):
            raise TypeError("oriented_length must be ScalarQuantity")
        if not isinstance(self.oriented_length.unit, PhysicalUnit):
            raise TypeError("oriented_length must use a physical unit")
        if not ModelSystemUnitConverter.compatible(
            self.oriented_length.unit,
            _METRE,
        ):
            raise ValueError("oriented_length must use units compatible with length")
        if self.oriented_length.magnitude == 0.0:
            raise ValueError("oriented_length must be nonzero")

    def _check_arg_unit_system(self) -> None:
        if type(self.unit_system) is not UnitSystem:
            raise TypeError("unit_system must be UnitSystem")
        if self.unit_system is UnitSystem.NONDIMENSIONAL:
            raise ValueError("unit_system must be physical")
        if self.oriented_length.unit != self.unit_system.length_unit:
            raise ValueError(
                "oriented_length must use the selected unit-system length unit"
            )

    @property
    def length(self) -> ScalarQuantity:
        """Return the positive primitive-cell length."""
        return ScalarQuantity(
            magnitude=abs(self.oriented_length.magnitude),
            unit=self.oriented_length.unit,
        )

    @property
    def reciprocal_basis(self) -> ScalarQuantity:
        r"""Return the oriented reciprocal basis $b=2\pi/a$."""
        return ScalarQuantity(
            magnitude=2.0 * np.pi / self.oriented_length.magnitude,
            unit=PhysicalUnit(f"1 / {self.oriented_length.unit.expression}"),
        )
