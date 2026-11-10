"""Tests for structured Slater--Koster spatial basis definitions."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.solidstate.tb.slater_koster.basis import (
    HomopolarSp3sStarSpatialBasis,
    RealCubicOrbitalComponentEnum,
    SlaterKosterBasisCenter,
    SlaterKosterSpatialBasisCoordinate,
    SlaterKosterSpatialBasisDefinition,
)
from projectkoios.physkit.units import UnitSystem, VectorQuantity


class TestHomopolarSp3sStarSpatialBasis:
    """Require exact ordered five-orbital blocks at declared centers."""

    def test__init__retains_ordered_centers_and_radial_channel(self) -> None:
        centers = (
            SlaterKosterBasisCenter(
                identifier="origin_center",
                fractional_coordinate=VectorQuantity(
                    magnitude=np.asarray((0.0, 0.0, 0.0)),
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                ),
            ),
            SlaterKosterBasisCenter(
                identifier="quarter_diagonal_center",
                fractional_coordinate=VectorQuantity(
                    magnitude=np.asarray((0.25, 0.25, 0.25)),
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                ),
            ),
        )
        components = (
            (RealCubicOrbitalComponentEnum.s, 0),
            (RealCubicOrbitalComponentEnum.p_x, 0),
            (RealCubicOrbitalComponentEnum.p_y, 0),
            (RealCubicOrbitalComponentEnum.p_z, 0),
            (RealCubicOrbitalComponentEnum.s, 1),
        )
        definition = SlaterKosterSpatialBasisDefinition(
            display_identifier="two-center-sp3s-star",
            centers=centers,
            coordinates=tuple(
                SlaterKosterSpatialBasisCoordinate(
                    center_identifier=center.identifier,
                    orbital_component=component,
                    radial_channel=radial_channel,
                )
                for center in centers
                for component, radial_channel in components
            ),
        )

        basis = HomopolarSp3sStarSpatialBasis(definition=definition)

        assert basis.dimension == 10
        assert basis.orbital_offset(center_identifier="origin_center") == 0
        assert basis.orbital_offset(center_identifier="quarter_diagonal_center") == 5
        np.testing.assert_array_equal(
            definition.orbital_centers_fractional.magnitude,
            np.vstack((np.zeros((5, 3)), np.full((5, 3), 0.25))),
        )
        assert definition.coordinates[4].radial_channel == 1
        assert (
            definition.coordinates[4].orbital_component
            is RealCubicOrbitalComponentEnum.s
        )

    def test__init__rejects_reordered_sp3s_star_block(self) -> None:
        center = SlaterKosterBasisCenter(
            identifier="origin_center",
            fractional_coordinate=VectorQuantity(
                magnitude=np.asarray((0.0, 0.0, 0.0)),
                unit=UnitSystem.NONDIMENSIONAL.length_unit,
            ),
        )
        coordinates = (
            SlaterKosterSpatialBasisCoordinate(
                center_identifier="origin_center",
                orbital_component=RealCubicOrbitalComponentEnum.p_x,
                radial_channel=0,
            ),
            SlaterKosterSpatialBasisCoordinate(
                center_identifier="origin_center",
                orbital_component=RealCubicOrbitalComponentEnum.s,
                radial_channel=0,
            ),
            SlaterKosterSpatialBasisCoordinate(
                center_identifier="origin_center",
                orbital_component=RealCubicOrbitalComponentEnum.p_y,
                radial_channel=0,
            ),
            SlaterKosterSpatialBasisCoordinate(
                center_identifier="origin_center",
                orbital_component=RealCubicOrbitalComponentEnum.p_z,
                radial_channel=0,
            ),
            SlaterKosterSpatialBasisCoordinate(
                center_identifier="origin_center",
                orbital_component=RealCubicOrbitalComponentEnum.s,
                radial_channel=1,
            ),
        )

        with pytest.raises(ValueError, match="exact ordered sp3s"):
            HomopolarSp3sStarSpatialBasis(
                definition=SlaterKosterSpatialBasisDefinition(
                    display_identifier="reordered",
                    centers=(center,),
                    coordinates=coordinates,
                )
            )
