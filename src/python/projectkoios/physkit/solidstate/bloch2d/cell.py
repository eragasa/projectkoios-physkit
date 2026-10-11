r"""Unit-aware two-dimensional direct cell for Bloch representations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.numerics.linear_algebra.volume import (
    columns_are_linearly_independent,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ModelSystemUnit,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    UnitSystem,
)

_NORMALIZED_AREA_ABS_TOLERANCE = 1.0e-8


def _squared_unit(unit: ModelSystemUnit) -> PhysicalUnit | Unitless:
    """Return the square of a supported direct-basis unit."""
    if isinstance(unit, Unitless):
        return Unitless()
    if isinstance(unit, PhysicalUnit):
        return PhysicalUnit(f"({unit.expression}) ** 2")
    raise TypeError("direct_basis unit must be PhysicalUnit or Unitless")


def _inverse_unit(unit: ModelSystemUnit) -> PhysicalUnit | Unitless:
    """Return the reciprocal of a supported direct-basis unit."""
    if isinstance(unit, Unitless):
        return Unitless()
    if isinstance(unit, PhysicalUnit):
        return PhysicalUnit(f"1 / {unit.expression}")
    raise TypeError("direct_basis unit must be PhysicalUnit or Unitless")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class BlochCell2D(DataObject):
    """Represent two direct-cell vectors stored as matrix columns.

    The cell may be represented in the metal unit system or as an explicitly
    nondimensional model. It is not required to be a primitive Bravais cell.
    """

    direct_basis: MatrixQuantity
    unit_system: UnitSystem

    def __post_init__(self) -> None:
        """Check every direct-cell argument."""
        if not isinstance(self.direct_basis, MatrixQuantity):
            raise TypeError("direct_basis must be MatrixQuantity")
        if type(self.unit_system) is not UnitSystem:
            raise TypeError("unit_system must be UnitSystem")
        if self.unit_system not in (
            UnitSystem.METAL,
            UnitSystem.NONDIMENSIONAL,
        ):
            raise ValueError("unit_system must be METAL or NONDIMENSIONAL")
        if self.direct_basis.unit != self.unit_system.length_unit:
            raise ValueError(
                "direct_basis must use the selected unit-system length unit"
            )
        if self.direct_basis.magnitude.shape != (2, 2):
            raise ValueError("direct_basis must have shape (2, 2)")
        if not columns_are_linearly_independent(
            self.direct_basis.magnitude,
            norm_vol_atol=_NORMALIZED_AREA_ABS_TOLERANCE,
        ):
            raise ValueError(
                "direct_basis columns must be numerically linearly independent"
            )

    @property
    def area(self) -> ScalarQuantity:
        """Return the positive direct-cell area."""
        return ScalarQuantity(
            magnitude=float(abs(np.linalg.det(self.direct_basis.magnitude))),
            unit=_squared_unit(self.direct_basis.unit),
        )

    @property
    def reciprocal_basis(self) -> MatrixQuantity:
        r"""Return $B=2\pi A^{-\mathsf T}$ in reciprocal basis units."""
        return MatrixQuantity(
            magnitude=np.asarray(
                2.0 * np.pi * np.linalg.inv(self.direct_basis.magnitude).T,
                dtype=np.float64,
            ),
            unit=_inverse_unit(self.direct_basis.unit),
        )

    @property
    def metric_tensor(self) -> MatrixQuantity:
        r"""Return $G=A^{\mathsf T}A$ in squared basis units."""
        basis = self.direct_basis.magnitude
        return MatrixQuantity(
            magnitude=np.asarray(basis.T @ basis, dtype=np.float64),
            unit=_squared_unit(self.direct_basis.unit),
        )
