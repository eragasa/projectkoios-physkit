r"""Construct spinless homopolar ``sp3s*`` two-center hopping matrices.

For source-row and target-column orbital order ``(s, p_x, p_y, p_z, s*)``, normalized
direction :math:`d=(l,m,n)`, and seven real energy integrals, the represented matrix is

.. math::

   T(d;V)=
   \begin{pmatrix}
   V_{ss\sigma} & V_{sp\sigma}d^{\mathsf T} & V_{ss^*\sigma}\\
   -V_{sp\sigma}d &
   V_{pp\pi}I_3+(V_{pp\sigma}-V_{pp\pi})dd^{\mathsf T} &
   -V_{s^*p\sigma}d\\
   V_{ss^*\sigma} & V_{s^*p\sigma}d^{\mathsf T} &
   V_{s^*s^*\sigma}
   \end{pmatrix}.

Equivalently, :math:`T=\sum_\alpha V_\alpha B_\alpha(d)`. The sign convention gives
:math:`T(-d;V)=T(d;V)^{\mathsf T}` for real homopolar integrals. This module owns only
the energy-valued two-center matrix. It owns no lattice, shell, material, onsite term,
overlap, spin, fitting policy, evidence, or scientific acceptance.
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
"""Fixed source-row and target-column order of the represented matrix."""


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarTwoCenterIntegrals:
    r"""Bind seven real physical energy integrals for one homopolar bond family.

    Parameters
    ----------
    ss_sigma, sp_sigma, pp_sigma, pp_pi, s_star_p_sigma, s_s_star_sigma,
    s_star_s_star_sigma
        Immutable scalar energy quantities in the stable channel order

        .. math::

           (ss\sigma,sp\sigma,pp\sigma,pp\pi,s^*p\sigma,ss^*\sigma,s^*s^*\sigma).

    unit_system
        Physical unit-system view whose exact energy unit is shared by all seven
        quantities. The default is `UnitSystem.ATOMIC`.

    Raises
    ------
    TypeError
        If the unit system, quantities, or quantity units lack exact nominal owners.
    ValueError
        If the selected system is nondimensional, a unit is not energy-compatible,
        units differ, or the common unit differs from the selected system energy unit.

    Notes
    -----
    An absent represented interaction is an explicit zero-valued energy quantity. It is
    never inferred from a neighbor-shell identifier. Valid storage does not establish
    fitting, transferability, convergence, uncertainty, or scientific acceptance.
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
        """Validate exact physical-energy ownership and one coherent unit."""
        if type(self.unit_system) is not UnitSystem:
            raise TypeError("unit_system must be UnitSystem")
        if self.unit_system is UnitSystem.NONDIMENSIONAL:
            raise ValueError("unit_system must provide physical energy units")

        values = self.ordered_quantities
        # Validate all channels, including explicit zeros, so a shell label can never
        # substitute for a dimensioned coefficient value.
        for name, value in zip(self.parameter_names, values, strict=True):
            self._validate_energy(value=value, name=name)

        # Exact equality retains one numerical scale. Compatible but different units
        # require an explicit conversion action before record construction.
        if any(value.unit != values[0].unit for value in values[1:]):
            raise ValueError("all two-center integrals must use the same energy unit")
        if values[0].unit != self.unit_system.energy_unit:
            raise ValueError("integrals must use the selected unit-system energy unit")

    @staticmethod
    def _validate_energy(value: ScalarQuantity, name: str) -> None:
        """Require one exact scalar quantity carrying a physical energy unit."""
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
        """Return stable channel names in exact integral order.

        Returns
        -------
        tuple[str, ...]
            Seven enum slug values aligned with `ordered_quantities`.
        """
        return tuple(channel.value for channel in HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER)

    @property
    def ordered_quantities(self) -> tuple[ScalarQuantity, ...]:
        """Return all seven integral quantities in stable channel order.

        Returns
        -------
        tuple[ScalarQuantity, ...]
            Immutable seven-element tuple used by channel-wise linear construction.
        """
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
        """Return the exact physical energy unit shared by every integral.

        Returns
        -------
        PhysicalUnit
            Common unit selected by ``unit_system``.

        Raises
        ------
        RuntimeError
            If an already validated immutable quantity is observed with changed unit
            ownership.
        """
        unit = self.ss_sigma.unit
        if type(unit) is not PhysicalUnit:
            raise RuntimeError("validated integral energy unit changed")
        return unit


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarTwoCenterHoppingDefinition:
    r"""Bind one source-to-target direction to one integral vector.

    Parameters
    ----------
    bond_direction
        Exact normalized Cartesian direction :math:`d` from source to target.
    integrals
        Exact seven-integral physical energy vector :math:`V`.

    Raises
    ------
    TypeError
        If either value lacks its exact nominal owner.

    Notes
    -----
    The pair :math:`(d,V)` fully defines the represented two-center matrix under this
    module's fixed orbital and sign convention. It does not define lattice support,
    center identities, onsite terms, or material interpretation.
    """

    bond_direction: SlaterKosterBondDirection
    integrals: HomopolarSp3sStarTwoCenterIntegrals

    def __post_init__(self) -> None:
        """Validate exact direction and integral-vector ownership."""
        if type(self.bond_direction) is not SlaterKosterBondDirection:
            raise TypeError("bond_direction must be SlaterKosterBondDirection")
        if type(self.integrals) is not HomopolarSp3sStarTwoCenterIntegrals:
            raise TypeError("integrals must be HomopolarSp3sStarTwoCenterIntegrals")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class HomopolarSp3sStarTwoCenterHopping(ResultsObject):
    r"""Correlate a two-center definition with its energy-valued matrix.

    Parameters
    ----------
    definition
        Exact direction and integral vector used to construct the matrix.
    matrix
        Immutable real ``(5, 5)`` `MatrixQuantity` in the integral energy unit.

    Raises
    ------
    TypeError
        If values lack their exact nominal owners.
    ValueError
        If the matrix shape or unit disagrees with the definition.

    Notes
    -----
    Rows are source orbitals and columns are target orbitals in
    `SP3S_STAR_ORBITAL_ORDER`. The result is one directed-bond matrix, not a complete
    Hamiltonian, fit result, or scientific validation record.
    """

    definition: HomopolarSp3sStarTwoCenterHoppingDefinition
    matrix: MatrixQuantity

    def __post_init__(self) -> None:
        """Validate exact definition, matrix shape, and energy unit correlation."""
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
        """Return the fixed source-row and target-column orbital order.

        Returns
        -------
        tuple[str, ...]
            ``("s", "p_x", "p_y", "p_z", "s_star")``.
        """
        return SP3S_STAR_ORBITAL_ORDER


class HomopolarSp3sStarTwoCenterHoppingConstructor:
    r"""Construct one real two-center matrix by channel-wise linear reduction.

    Notes
    -----
    The constructor delegates angular placement to
    `HomopolarSp3sStarTwoCenterChannelBasisConstructor` and evaluates

    .. math::

       T(d;V)=\sum_{\alpha=1}^{7}V_\alpha B_\alpha(d)

    in stable channel order. It performs no shell inference, fitting, regularization,
    onsite construction, lattice assembly, or scientific comparison.
    """

    __slots__ = ()

    def execute(
        self,
        *,
        definition: HomopolarSp3sStarTwoCenterHoppingDefinition,
    ) -> HomopolarSp3sStarTwoCenterHopping:
        r"""Construct the immutable energy-valued hopping matrix.

        Parameters
        ----------
        definition
            Exact direction and seven-integral vector. The argument is keyword-only.

        Returns
        -------
        HomopolarSp3sStarTwoCenterHopping
            Correlated definition and immutable real ``(5, 5)`` energy matrix.

        Raises
        ------
        TypeError
            If ``definition`` lacks its exact nominal owner.
        OverflowError
            If finite input magnitudes produce any nonfinite matrix entry during
            binary64 accumulation. Values are not clipped or rescaled.

        Notes
        -----
        The represented sign convention gives

        .. math::

           T(-d;V)=T(d;V)^{\mathsf T}

        for real homopolar integrals.
        """
        if type(definition) is not HomopolarSp3sStarTwoCenterHoppingDefinition:
            raise TypeError(
                "definition must be HomopolarSp3sStarTwoCenterHoppingDefinition"
            )
        integrals = definition.integrals

        # Angular matrices are nondimensional and ordered exactly like the integral
        # tuple, so the reduction requires no string matching or coefficient inference.
        channel_basis = HomopolarSp3sStarTwoCenterChannelBasisConstructor().execute(
            bond_direction=definition.bond_direction
        )
        matrix: Float64Matrix = np.zeros((5, 5), dtype=np.float64)

        # Suppress intermediate NumPy warnings so overflow is reported through the
        # explicit public exception below rather than as warning-dependent behavior.
        with np.errstate(over="ignore", invalid="ignore"):
            for coefficient, basis_matrix in zip(
                integrals.ordered_quantities,
                channel_basis.matrices,
                strict=True,
            ):
                matrix += coefficient.magnitude * basis_matrix.magnitude

        # Finite quantities can still overflow a finite-precision linear combination;
        # fail closed instead of returning clipped or nonfinite represented physics.
        if not np.all(np.isfinite(matrix)):
            raise OverflowError("two-center hopping arithmetic overflowed")
        return HomopolarSp3sStarTwoCenterHopping(
            definition=definition,
            matrix=MatrixQuantity(
                magnitude=matrix,
                unit=integrals.energy_unit,
            ),
        )
