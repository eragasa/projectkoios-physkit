"""Public package façade checks for Slater--Koster mechanics."""

from __future__ import annotations

import projectkoios.physkit.solidstate.tb.slater_koster as slater_koster
from projectkoios.physkit.solidstate.tb.slater_koster.sp3s_star import (
    HomopolarSp3sStarTwoCenterHoppingConstructor,
)


class TestSlaterKosterPublicApi:
    """Keep the coherent specialization surface deliberate and explicit."""

    def test__exports_two_center_constructor_from_defining_module(self) -> None:
        assert (
            slater_koster.HomopolarSp3sStarTwoCenterHoppingConstructor
            is HomopolarSp3sStarTwoCenterHoppingConstructor
        )
        assert slater_koster.__all__ == [
            "HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER",
            "SP3S_STAR_ORBITAL_ORDER",
            "HomopolarSp3sStarBondGeometry",
            "HomopolarSp3sStarIntegralChannelEnum",
            "HomopolarSp3sStarParameterGroup",
            "HomopolarSp3sStarParameterizedBond",
            "HomopolarSp3sStarSpatialBasis",
            "HomopolarSp3sStarTwoCenterChannelBasis",
            "HomopolarSp3sStarTwoCenterChannelBasisConstructor",
            "HomopolarSp3sStarTwoCenterHopping",
            "HomopolarSp3sStarTwoCenterHoppingConstructor",
            "HomopolarSp3sStarTwoCenterHoppingDefinition",
            "HomopolarSp3sStarTwoCenterIntegrals",
            "RealCubicOrbitalComponentEnum",
            "SlaterKosterBasisCenter",
            "SlaterKosterBondDirection",
            "SlaterKosterDirectedBond",
            "SlaterKosterSpatialBasisCoordinate",
            "SlaterKosterSpatialBasisDefinition",
        ]
