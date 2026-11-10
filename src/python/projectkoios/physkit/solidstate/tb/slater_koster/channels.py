r"""Construct the linear angular-channel basis for homopolar ``sp3s*`` hopping.

For normalized source-to-target direction :math:`d=(l,m,n)` and the fixed orbital order
``(s, p_x, p_y, p_z, s*)``, the energy-valued two-center matrix is linear in seven
integrals:

.. math::

   T(d;V)=\sum_{\alpha=1}^{7}V_\alpha B_\alpha(d).

This module owns the stable channel enumeration and the seven nondimensional real
coefficient matrices :math:`B_\alpha`. It owns no energy values, material semantics,
fitting, lattice support, onsite terms, or scientific acceptance.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from projectkoios.numerics.numpy.typing import Float64Matrix
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import MatrixQuantity, UnitSystem

from .bond import SlaterKosterBondDirection


class HomopolarSp3sStarIntegralChannelEnum(StrEnum):
    r"""Enumerate the stable seven-channel homopolar ``sp3s*`` vocabulary.

    Notes
    -----
    The enum order is observable and shared by integral records, channel matrices, and
    the direct hopping constructor:

    .. math::

       (ss\sigma,sp\sigma,pp\sigma,pp\pi,s^*p\sigma,ss^*\sigma,s^*s^*\sigma).

    This closed vocabulary does not claim support for heteropolar, overlap,
    spin-dependent, ``d``-orbital, or ``f``-orbital interaction tables.
    """

    ss_sigma = "ss_sigma"
    sp_sigma = "sp_sigma"
    pp_sigma = "pp_sigma"
    pp_pi = "pp_pi"
    s_star_p_sigma = "s_star_p_sigma"
    s_s_star_sigma = "s_s_star_sigma"
    s_star_s_star_sigma = "s_star_s_star_sigma"


HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER = tuple(HomopolarSp3sStarIntegralChannelEnum)
"""Stable channel order used by integral, matrix, and later linear records."""


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class HomopolarSp3sStarTwoCenterChannelBasis(ResultsObject):
    r"""Correlate one bond direction with seven angular coefficient matrices.

    Parameters
    ----------
    bond_direction
        Exact normalized source-to-target Cartesian direction.
    matrices
        Seven-element tuple of immutable nondimensional real ``(5, 5)`` matrix
        quantities in `HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER`.

    Raises
    ------
    TypeError
        If values lack their exact nominal owners.
    ValueError
        If the tuple length, units, or matrix shapes differ from the seven-channel
        contract.

    Notes
    -----
    The result represents

    .. math::

       \mathcal B_{\mathrm{channel}}(d)
       =(B_1(d),\ldots,B_7(d)).

    It contains no energy coefficients and does not establish parameter identifiability
    or material sufficiency.
    """

    bond_direction: SlaterKosterBondDirection
    matrices: tuple[MatrixQuantity, ...]

    def __post_init__(self) -> None:
        """Validate direction ownership and exact seven-matrix structure."""
        if type(self.bond_direction) is not SlaterKosterBondDirection:
            raise TypeError("bond_direction must be SlaterKosterBondDirection")
        if type(self.matrices) is not tuple:
            raise TypeError("matrices must be a tuple")
        if len(self.matrices) != len(HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER):
            raise ValueError("matrices must contain exactly seven channel matrices")
        if any(type(matrix) is not MatrixQuantity for matrix in self.matrices):
            raise TypeError("matrices must contain MatrixQuantity values")

        # Angular coefficients are pure numbers. Retaining the nondimensional unit
        # makes that convention explicit without introducing a parallel unit API.
        if any(
            matrix.unit != UnitSystem.NONDIMENSIONAL.length_unit
            for matrix in self.matrices
        ):
            raise ValueError("channel matrices must be nondimensional")
        if any(matrix.magnitude.shape != (5, 5) for matrix in self.matrices):
            raise ValueError("channel matrices must have shape (5, 5)")

    def matrix(
        self,
        *,
        channel: HomopolarSp3sStarIntegralChannelEnum,
    ) -> MatrixQuantity:
        """Return the exact coefficient matrix for one channel.

        Parameters
        ----------
        channel
            Exact channel enum member. The argument is keyword-only.

        Returns
        -------
        MatrixQuantity
            Immutable nondimensional ``(5, 5)`` coefficient matrix.

        Raises
        ------
        TypeError
            If ``channel`` is not exactly `HomopolarSp3sStarIntegralChannelEnum`.
        """
        if type(channel) is not HomopolarSp3sStarIntegralChannelEnum:
            raise TypeError("channel must be HomopolarSp3sStarIntegralChannelEnum")

        # Enum order is the single indexing convention shared across all seven-channel
        # records; no string lookup or inferred shell semantics is used.
        return self.matrices[HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER.index(channel)]


class HomopolarSp3sStarTwoCenterChannelBasisConstructor:
    r"""Construct seven nondimensional angular matrices from one direction.

    Notes
    -----
    With :math:`B_\alpha=\partial T/\partial V_\alpha`, every matrix is independent of
    energy values. The implementation uses the fixed orbital order
    ``(s, p_x, p_y, p_z, s*)`` and the sign convention

    .. math::

       (B_{sp\sigma})_{s,p_i}=d_i,
       \qquad
       (B_{sp\sigma})_{p_i,s}=-d_i.
    """

    __slots__ = ()

    def execute(
        self,
        *,
        bond_direction: SlaterKosterBondDirection,
    ) -> HomopolarSp3sStarTwoCenterChannelBasis:
        r"""Construct the exact ordered seven-matrix channel basis.

        Parameters
        ----------
        bond_direction
            Normalized source-to-target direction. The argument is keyword-only.

        Returns
        -------
        HomopolarSp3sStarTwoCenterChannelBasis
            Correlated direction and immutable matrix tuple.

        Raises
        ------
        TypeError
            If ``bond_direction`` lacks its exact nominal owner.

        Notes
        -----
        Construction uses binary64 arrays. Each branch writes only entries owned by the
        corresponding channel; all other entries remain exact zeros.
        """
        if type(bond_direction) is not SlaterKosterBondDirection:
            raise TypeError("bond_direction must be SlaterKosterBondDirection")
        direction = bond_direction.components.magnitude
        matrices: list[MatrixQuantity] = []

        # Iterating the enum-derived tuple guarantees that output matrix order matches
        # integral order without a second mapping table.
        for channel in HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER:
            matrix: Float64Matrix = np.zeros((5, 5), dtype=np.float64)

            if channel is HomopolarSp3sStarIntegralChannelEnum.ss_sigma:
                # Scalar s-to-s channel.
                matrix[0, 0] = 1.0
            elif channel is HomopolarSp3sStarIntegralChannelEnum.sp_sigma:
                # Antisymmetric source/target s-p sign convention.
                matrix[0, 1:4] = direction
                matrix[1:4, 0] = -direction
            elif channel is HomopolarSp3sStarIntegralChannelEnum.pp_sigma:
                # Longitudinal projector d d^T in the p subspace.
                matrix[1:4, 1:4] = np.outer(direction, direction)
            elif channel is HomopolarSp3sStarIntegralChannelEnum.pp_pi:
                # Transverse projector I - d d^T in the p subspace.
                matrix[1:4, 1:4] = np.eye(3) - np.outer(direction, direction)
            elif channel is HomopolarSp3sStarIntegralChannelEnum.s_star_p_sigma:
                # The s-star coordinate is index four, with the same orientation
                # convention as the primary s-p channel.
                matrix[4, 1:4] = direction
                matrix[1:4, 4] = -direction
            elif channel is HomopolarSp3sStarIntegralChannelEnum.s_s_star_sigma:
                # Symmetric coupling between radial s channels zero and one.
                matrix[0, 4] = 1.0
                matrix[4, 0] = 1.0
            elif channel is HomopolarSp3sStarIntegralChannelEnum.s_star_s_star_sigma:
                # Scalar s-star-to-s-star channel.
                matrix[4, 4] = 1.0

            # Quantity construction freezes each completed angular matrix before it is
            # exposed through the correlated result.
            matrices.append(
                MatrixQuantity(
                    magnitude=matrix,
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                )
            )
        return HomopolarSp3sStarTwoCenterChannelBasis(
            bond_direction=bond_direction,
            matrices=tuple(matrices),
        )
