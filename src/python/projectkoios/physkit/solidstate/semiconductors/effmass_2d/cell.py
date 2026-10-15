r"""Unit-aware two-dimensional periodic primitive-cell geometry."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    UnitSystem,
)

_METRE = PhysicalUnit(expression="meter")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicPrimitiveCell2D(DataObject):
    """Represent two physical primitive vectors stored as matrix columns."""

    primitive_basis: MatrixQuantity
    unit_system: UnitSystem

    def __post_init__(self) -> None:
        """Check every primitive-cell argument."""
        self._check_arg_primitive_basis()
        self._check_arg_unit_system()

    def _check_arg_primitive_basis(self) -> None:
        """Require an invertible two-dimensional physical length basis."""
        if not isinstance(self.primitive_basis, MatrixQuantity):
            raise TypeError("primitive_basis must be MatrixQuantity")
        if not isinstance(self.primitive_basis.unit, PhysicalUnit):
            raise TypeError("primitive_basis must use a physical length unit")
        if not MODEL_SYSTEM_UNIT_CONVERTER.compatible(
            self.primitive_basis.unit,
            _METRE,
        ):
            raise ValueError("primitive_basis must use units compatible with length")
        if self.primitive_basis.magnitude.shape != (2, 2):
            raise ValueError("primitive_basis must have shape (2, 2)")
        if np.linalg.det(self.primitive_basis.magnitude) == 0.0:
            raise ValueError("primitive_basis columns must be linearly independent")

    def _check_arg_unit_system(self) -> None:
        """Require a physical coherent unit system matching the basis unit."""
        if type(self.unit_system) is not UnitSystem:
            raise TypeError("unit_system must be UnitSystem")
        if self.unit_system is UnitSystem.NONDIMENSIONAL:
            raise ValueError("unit_system must be physical")
        if self.primitive_basis.unit != self.unit_system.length_unit:
            raise ValueError(
                "primitive_basis must use the selected unit-system length unit"
            )

    @property
    def area(self) -> ScalarQuantity:
        """Return the positive primitive-cell area in the squared input unit."""
        return ScalarQuantity(
            magnitude=float(abs(np.linalg.det(self.primitive_basis.magnitude))),
            unit=PhysicalUnit(
                expression=f"({self.primitive_basis.unit.expression}) ** 2"
            ),
        )

    @property
    def reciprocal_basis(self) -> MatrixQuantity:
        r"""Return $B=2\pi A^{-\mathsf T}$ in inverse input-length units."""
        return MatrixQuantity(
            magnitude=np.asarray(
                2.0 * np.pi * np.linalg.inv(self.primitive_basis.magnitude).T,
                dtype=np.float64,
            ),
            unit=PhysicalUnit(expression=f"1 / {self.primitive_basis.unit.expression}"),
        )

    @property
    def metric_tensor(self) -> MatrixQuantity:
        r"""Return $G=A^{\mathsf T}A$ in squared input-length units."""
        basis = self.primitive_basis.magnitude
        return MatrixQuantity(
            magnitude=np.asarray(basis.T @ basis, dtype=np.float64),
            unit=PhysicalUnit(
                expression=f"({self.primitive_basis.unit.expression}) ** 2"
            ),
        )
