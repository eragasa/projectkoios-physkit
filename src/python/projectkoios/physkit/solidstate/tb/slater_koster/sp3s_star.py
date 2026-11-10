r"""Spinless homopolar ``sp3s*`` two-center Slater--Koster construction.

This module owns only the bond-direction map from seven explicit
Slater--Koster energy integrals to one real five-by-five hopping matrix. It owns no
lattice, neighbor shell, Bloch phase, material selection, fitting policy, evidence, or
scientific acceptance.

For the unit direction :math:`d=(l,m,n)` from the row/source center to the
column/target center, the fixed orbital order is ``(s, p_x, p_y, p_z, s*)`` and the
sign convention is

.. math::

   \langle s|H|p_i\rangle=d_i V_{sp\sigma},\qquad
   \langle p_i|H|s\rangle=-d_i V_{sp\sigma}.

The same convention applies to the ``s*``--``p`` channel. Consequently
:math:`T(-d)=T(d)^\mathsf{T}` for real homopolar integrals.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.numerics.numpy.typing import Float64Matrix
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    PintUnitConverter,
    ScalarQuantity,
    UnitSystem,
)

from .bond import SlaterKosterBondDirection
from .channels import (
    HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER,
    HomopolarSp3sStarTwoCenterChannelBasisConstructor,
)

SP3S_STAR_ORBITAL_ORDER = ("s", "p_x", "p_y", "p_z", "s_star")
"""Fixed row and column order of the represented two-center matrix."""


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarTwoCenterIntegrals:
    r"""Bind seven real energy integrals for one homopolar ``sp3s*`` bond.

    All values must carry the same exact PhysKit unit. Zero-valued absent channels are
    explicit quantities rather than inferred from a neighbor-shell label.
    """

    ss_sigma: ScalarQuantity
    sp_sigma: ScalarQuantity
    pp_sigma: ScalarQuantity
    pp_pi: ScalarQuantity
    s_star_p_sigma: ScalarQuantity
    s_s_star_sigma: ScalarQuantity
    s_star_s_star_sigma: ScalarQuantity
    unit_system: UnitSystem = UnitSystem.ATOMIC

    def __post_init__(self) -> None:
        if type(self.unit_system) is not UnitSystem:
            raise TypeError("unit_system must be UnitSystem")
        if self.unit_system is UnitSystem.NONDIMENSIONAL:
            raise ValueError("unit_system must provide physical energy units")
        values = self.ordered_quantities
        for name, value in zip(self.parameter_names, values, strict=True):
            self._validate_energy(value=value, name=name)
        if any(value.unit != values[0].unit for value in values[1:]):
            raise ValueError("all two-center integrals must use the same energy unit")
        if values[0].unit != self.unit_system.energy_unit:
            raise ValueError("integrals must use the selected unit-system energy unit")

    @staticmethod
    def _validate_energy(value: ScalarQuantity, name: str) -> None:
        """Require one scalar carrying a physical energy unit."""
        if type(value) is not ScalarQuantity:
            raise TypeError(f"{name} must be ScalarQuantity")
        if type(value.unit) is not PhysicalUnit:
            raise TypeError(f"{name}.unit must be PhysicalUnit")
        if not PintUnitConverter().compatible(
            left=value.unit,
            right=PhysicalUnit(expression="joule"),
        ):
            raise ValueError(f"{name}.unit must represent energy")

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Return the stable seven-channel order."""
        return tuple(channel.value for channel in HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER)

    @property
    def ordered_quantities(self) -> tuple[ScalarQuantity, ...]:
        """Return the integral quantities in stable channel order."""
        return (
            self.ss_sigma,
            self.sp_sigma,
            self.pp_sigma,
            self.pp_pi,
            self.s_star_p_sigma,
            self.s_s_star_sigma,
            self.s_star_s_star_sigma,
        )

    @property
    def energy_unit(self) -> PhysicalUnit:
        """Return the exact physical energy unit shared by every integral."""
        unit = self.ss_sigma.unit
        if type(unit) is not PhysicalUnit:
            raise RuntimeError("validated integral energy unit changed")
        return unit


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarTwoCenterHoppingDefinition:
    """Bind a unit direction from source to target and its seven bond integrals."""

    bond_direction: SlaterKosterBondDirection
    integrals: HomopolarSp3sStarTwoCenterIntegrals

    def __post_init__(self) -> None:
        if type(self.bond_direction) is not SlaterKosterBondDirection:
            raise TypeError("bond_direction must be SlaterKosterBondDirection")
        if type(self.integrals) is not HomopolarSp3sStarTwoCenterIntegrals:
            raise TypeError("integrals must be HomopolarSp3sStarTwoCenterIntegrals")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class HomopolarSp3sStarTwoCenterHopping(ResultsObject):
    """Retain one definition and its real five-orbital hopping matrix."""

    definition: HomopolarSp3sStarTwoCenterHoppingDefinition
    matrix: MatrixQuantity

    def __post_init__(self) -> None:
        if type(self.definition) is not HomopolarSp3sStarTwoCenterHoppingDefinition:
            raise TypeError(
                "definition must be HomopolarSp3sStarTwoCenterHoppingDefinition"
            )
        if type(self.matrix) is not MatrixQuantity:
            raise TypeError("matrix must be MatrixQuantity")
        if self.matrix.magnitude.shape != (5, 5):
            raise ValueError("matrix must have shape (5, 5)")
        if self.matrix.unit != self.definition.integrals.energy_unit:
            raise ValueError("matrix unit must match the integral energy unit")

    @property
    def orbital_order(self) -> tuple[str, ...]:
        """Return the fixed source-row and target-column orbital order."""
        return SP3S_STAR_ORBITAL_ORDER


class HomopolarSp3sStarTwoCenterHoppingConstructor:
    """Construct the real two-center matrix from a bond direction and integrals."""

    __slots__ = ()

    def execute(
        self,
        *,
        definition: HomopolarSp3sStarTwoCenterHoppingDefinition,
    ) -> HomopolarSp3sStarTwoCenterHopping:
        """Return the represented hopping matrix without shell or lattice inference.

        Arithmetic uses binary64. Finite inputs that overflow matrix construction raise
        :class:`OverflowError`; values are not clipped or coerced back to finite data.
        """
        if type(definition) is not HomopolarSp3sStarTwoCenterHoppingDefinition:
            raise TypeError(
                "definition must be HomopolarSp3sStarTwoCenterHoppingDefinition"
            )
        integrals = definition.integrals
        channel_basis = HomopolarSp3sStarTwoCenterChannelBasisConstructor().execute(
            bond_direction=definition.bond_direction
        )
        matrix: Float64Matrix = np.zeros((5, 5), dtype=np.float64)
        with np.errstate(over="ignore", invalid="ignore"):
            for coefficient, basis_matrix in zip(
                integrals.ordered_quantities,
                channel_basis.matrices,
                strict=True,
            ):
                matrix += coefficient.magnitude * basis_matrix.magnitude
        if not np.all(np.isfinite(matrix)):
            raise OverflowError("two-center hopping arithmetic overflowed")
        return HomopolarSp3sStarTwoCenterHopping(
            definition=definition,
            matrix=MatrixQuantity(
                magnitude=matrix,
                unit=integrals.energy_unit,
            ),
        )
