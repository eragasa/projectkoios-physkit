"""Endpoint-excluded grids on one physical periodic cell."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    VectorQuantity,
)

_METRE = PhysicalUnit("meter")


@dataclass(frozen=True, slots=True, kw_only=True)
class PeriodicFiniteDifferenceGrid1D(DataObject):
    """Define an endpoint-excluded uniform grid on one oriented cell."""

    point_count: int
    oriented_cell_length: ScalarQuantity

    def __post_init__(self) -> None:
        """Check every grid argument."""
        self._check_arg_point_count()
        self._check_arg_oriented_cell_length()

    def _check_arg_point_count(self) -> None:
        if type(self.point_count) is not int:
            raise TypeError("point_count must be a built-in integer")
        if self.point_count < 3:
            raise ValueError("point_count must be at least three")

    def _check_arg_oriented_cell_length(self) -> None:
        if not isinstance(self.oriented_cell_length, ScalarQuantity):
            raise TypeError("oriented_cell_length must be ScalarQuantity")
        if not isinstance(self.oriented_cell_length.unit, PhysicalUnit):
            raise TypeError("oriented_cell_length must use a physical unit")
        if not ModelSystemUnitConverter.compatible(
            self.oriented_cell_length.unit,
            _METRE,
        ):
            raise ValueError(
                "oriented_cell_length must use units compatible with length"
            )
        if self.oriented_cell_length.magnitude == 0.0:
            raise ValueError("oriented_cell_length must be nonzero")

    @property
    def spacing(self) -> ScalarQuantity:
        """Return the positive endpoint-excluded grid spacing."""
        return ScalarQuantity(
            magnitude=abs(self.oriented_cell_length.magnitude) / self.point_count,
            unit=self.oriented_cell_length.unit,
        )

    @property
    def positions(self) -> VectorQuantity:
        """Return oriented cell positions without the repeated endpoint."""
        step = self.oriented_cell_length.magnitude / self.point_count
        return VectorQuantity(
            magnitude=np.arange(self.point_count, dtype=np.float64) * step,
            unit=self.oriented_cell_length.unit,
        )
