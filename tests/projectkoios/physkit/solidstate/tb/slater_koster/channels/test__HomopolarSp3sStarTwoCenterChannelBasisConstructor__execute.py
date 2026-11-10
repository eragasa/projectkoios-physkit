"""Tests for the homopolar ``sp3s*`` angular channel basis."""

from __future__ import annotations

import numpy as np

from projectkoios.physkit.solidstate.tb.slater_koster.bond import (
    SlaterKosterBondDirection,
)
from projectkoios.physkit.solidstate.tb.slater_koster.channels import (
    HomopolarSp3sStarIntegralChannelEnum,
    HomopolarSp3sStarTwoCenterChannelBasisConstructor,
)
from projectkoios.physkit.units import UnitSystem, VectorQuantity


class TestHomopolarSp3sStarTwoCenterChannelBasisConstructor:
    """Verify every nondimensional channel matrix independently."""

    def test__execute__constructs_all_seven_channel_matrices(self) -> None:
        direction = np.asarray((1.0, -2.0, 2.0), dtype=np.float64) / 3.0
        bond_direction = SlaterKosterBondDirection(
            components=VectorQuantity(
                magnitude=direction,
                unit=UnitSystem.NONDIMENSIONAL.length_unit,
            )
        )

        result = HomopolarSp3sStarTwoCenterChannelBasisConstructor().execute(
            bond_direction=bond_direction
        )

        ss = np.zeros((5, 5))
        ss[0, 0] = 1.0
        sp = np.zeros((5, 5))
        sp[0, 1:4] = direction
        sp[1:4, 0] = -direction
        pp_sigma = np.zeros((5, 5))
        pp_sigma[1:4, 1:4] = np.outer(direction, direction)
        pp_pi = np.zeros((5, 5))
        pp_pi[1:4, 1:4] = np.eye(3) - np.outer(direction, direction)
        s_star_p = np.zeros((5, 5))
        s_star_p[4, 1:4] = direction
        s_star_p[1:4, 4] = -direction
        s_s_star = np.zeros((5, 5))
        s_s_star[0, 4] = 1.0
        s_s_star[4, 0] = 1.0
        s_star_s_star = np.zeros((5, 5))
        s_star_s_star[4, 4] = 1.0

        expected = {
            HomopolarSp3sStarIntegralChannelEnum.ss_sigma: ss,
            HomopolarSp3sStarIntegralChannelEnum.sp_sigma: sp,
            HomopolarSp3sStarIntegralChannelEnum.pp_sigma: pp_sigma,
            HomopolarSp3sStarIntegralChannelEnum.pp_pi: pp_pi,
            HomopolarSp3sStarIntegralChannelEnum.s_star_p_sigma: s_star_p,
            HomopolarSp3sStarIntegralChannelEnum.s_s_star_sigma: s_s_star,
            HomopolarSp3sStarIntegralChannelEnum.s_star_s_star_sigma: s_star_s_star,
        }
        for channel, expected_matrix in expected.items():
            matrix = result.matrix(channel=channel)
            np.testing.assert_array_equal(matrix.magnitude, expected_matrix)
            assert matrix.unit == UnitSystem.NONDIMENSIONAL.length_unit
            assert matrix.magnitude.flags.writeable is False
