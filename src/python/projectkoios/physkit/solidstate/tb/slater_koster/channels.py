"""Linear angular-channel basis for homopolar ``sp3s*`` hopping."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from projectkoios.numerics.numpy.typing import Float64Matrix
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import MatrixQuantity, UnitSystem

from .bond import SlaterKosterBondDirection


class HomopolarSp3sStarIntegralChannelEnum(StrEnum):
    """Stable seven-channel order for the current homopolar specialization."""

    ss_sigma = "ss_sigma"
    sp_sigma = "sp_sigma"
    pp_sigma = "pp_sigma"
    pp_pi = "pp_pi"
    s_star_p_sigma = "s_star_p_sigma"
    s_s_star_sigma = "s_s_star_sigma"
    s_star_s_star_sigma = "s_star_s_star_sigma"


HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER = tuple(HomopolarSp3sStarIntegralChannelEnum)
"""Stable channel order used by integral and linear-parameter records."""


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class HomopolarSp3sStarTwoCenterChannelBasis(ResultsObject):
    """Retain seven ordered nondimensional angular coefficient matrices."""

    bond_direction: SlaterKosterBondDirection
    matrices: tuple[MatrixQuantity, ...]

    def __post_init__(self) -> None:
        if type(self.bond_direction) is not SlaterKosterBondDirection:
            raise TypeError("bond_direction must be SlaterKosterBondDirection")
        if type(self.matrices) is not tuple:
            raise TypeError("matrices must be a tuple")
        if len(self.matrices) != len(HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER):
            raise ValueError("matrices must contain exactly seven channel matrices")
        if any(type(matrix) is not MatrixQuantity for matrix in self.matrices):
            raise TypeError("matrices must contain MatrixQuantity values")
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
        """Return the exact coefficient matrix for one explicit channel."""
        if type(channel) is not HomopolarSp3sStarIntegralChannelEnum:
            raise TypeError("channel must be HomopolarSp3sStarIntegralChannelEnum")
        return self.matrices[HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER.index(channel)]


class HomopolarSp3sStarTwoCenterChannelBasisConstructor:
    """Construct seven nondimensional matrices from one bond direction."""

    __slots__ = ()

    def execute(
        self,
        *,
        bond_direction: SlaterKosterBondDirection,
    ) -> HomopolarSp3sStarTwoCenterChannelBasis:
        """Return the exact linear channel basis without energy coefficients."""
        if type(bond_direction) is not SlaterKosterBondDirection:
            raise TypeError("bond_direction must be SlaterKosterBondDirection")
        direction = bond_direction.components.magnitude
        matrices: list[MatrixQuantity] = []
        for channel in HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER:
            matrix: Float64Matrix = np.zeros((5, 5), dtype=np.float64)
            if channel is HomopolarSp3sStarIntegralChannelEnum.ss_sigma:
                matrix[0, 0] = 1.0
            elif channel is HomopolarSp3sStarIntegralChannelEnum.sp_sigma:
                matrix[0, 1:4] = direction
                matrix[1:4, 0] = -direction
            elif channel is HomopolarSp3sStarIntegralChannelEnum.pp_sigma:
                matrix[1:4, 1:4] = np.outer(direction, direction)
            elif channel is HomopolarSp3sStarIntegralChannelEnum.pp_pi:
                matrix[1:4, 1:4] = np.eye(3) - np.outer(direction, direction)
            elif channel is HomopolarSp3sStarIntegralChannelEnum.s_star_p_sigma:
                matrix[4, 1:4] = direction
                matrix[1:4, 4] = -direction
            elif channel is HomopolarSp3sStarIntegralChannelEnum.s_s_star_sigma:
                matrix[0, 4] = 1.0
                matrix[4, 0] = 1.0
            elif channel is HomopolarSp3sStarIntegralChannelEnum.s_star_s_star_sigma:
                matrix[4, 4] = 1.0
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
