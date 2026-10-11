r"""Unit-aware homogeneous-Dirichlet particle in a parallelogram."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.numerics.linear_algebra.volume import (
    columns_are_linearly_independent,
)
from projectkoios.physkit.periodic.lattice.lattice2d import DirectLattice2D
from projectkoios.physkit.periodic.lattice.metric import DirectLatticeMetric
from projectkoios.physkit.qm.models.model2d import BaseQuantumModel2D
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    UnitSystem,
)

_NORMALIZED_AREA_ABS_TOLERANCE = 1.0e-8


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class ObliquePiab2D(BaseQuantumModel2D):
    r"""Define a free particle confined to a direct-basis parallelogram.

    Cartesian position is :math:`\mathbf r=A\mathbf s`, where the columns of
    ``direct_basis`` are the two parallelogram edges and
    :math:`\mathbf s\in[0,1]^2`. The wavefunction obeys homogeneous Dirichlet
    conditions on all four edges. This model is bounded and is not a periodic
    primitive-cell or Bloch model.

    Parameters
    ----------
    direct_basis
        Finite linearly independent edge vectors stored as columns of a
        ``(2, 2)`` matrix. Its unit must equal ``unit_system.length_unit``.
    mass
        Finite positive particle mass in ``unit_system.mass_unit``.
    unit_system
        Coherent numerical unit system for length, mass, action, and energy.

    Raises
    ------
    TypeError
        If an argument has the wrong nominal runtime type.
    ValueError
        If units, shape, rank, finiteness, or mass positivity are invalid.
    """

    direct_basis: MatrixQuantity
    mass: ScalarQuantity
    unit_system: UnitSystem

    def __post_init__(self) -> None:
        """Check geometry, mass, and coherent-unit declarations."""
        if not isinstance(self.direct_basis, MatrixQuantity):
            raise TypeError("direct_basis must be MatrixQuantity")
        if not isinstance(self.mass, ScalarQuantity):
            raise TypeError("mass must be ScalarQuantity")
        if type(self.unit_system) is not UnitSystem:
            raise TypeError("unit_system must be UnitSystem")
        if self.direct_basis.unit != self.unit_system.length_unit:
            raise ValueError("direct_basis must use the selected length unit")
        if self.direct_basis.magnitude.shape != (2, 2):
            raise ValueError("direct_basis must have shape (2, 2)")
        if not columns_are_linearly_independent(
            self.direct_basis.magnitude,
            norm_vol_atol=_NORMALIZED_AREA_ABS_TOLERANCE,
        ):
            raise ValueError("direct_basis columns must be linearly independent")
        if self.mass.unit != self.unit_system.mass_unit:
            raise ValueError("mass must use the selected mass unit")
        if self.mass.magnitude <= 0.0:
            raise ValueError("mass must be positive")

    @property
    def direct_lattice(self) -> DirectLattice2D:
        """Return the unit-free direct geometry used for metric operations."""
        basis = self.direct_basis.magnitude
        return DirectLattice2D(a1=basis[:, 0], a2=basis[:, 1])

    @property
    def metric(self) -> DirectLatticeMetric:
        """Return the constant direct metric and its inverse."""
        return DirectLatticeMetric(self.direct_lattice)

    @property
    def area(self) -> ScalarQuantity:
        """Return the positive physical parallelogram area."""
        if self.unit_system is UnitSystem.NONDIMENSIONAL:
            unit = Unitless()
        else:
            unit = PhysicalUnit(f"({self.direct_basis.unit.expression}) ** 2")
        return ScalarQuantity(
            magnitude=float(abs(np.linalg.det(self.direct_basis.magnitude))),
            unit=unit,
        )
