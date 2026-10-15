"""Argument checks for empirical carrier mobility."""

import pytest

from projectkoios.physkit.constants import SI
from projectkoios.physkit.solidstate.semiconductors.mobility import (
    CarrierMobilityModel,
)


def test__CarrierMobilityModel__init__rejects_nonpositive_reference_data() -> None:
    with pytest.raises(ValueError, match="lattice_mobility_reference"):
        CarrierMobilityModel(
            effective_mass=0.26 * SI.me0,
            lattice_mobility_reference=0.0,
            impurity_mobility_reference=1_500.0,
            impurity_concentration_reference=1.0e16,
        )


def test__CarrierMobilityModel__lattice_mobility__rejects_nonpositive_temperature() -> (
    None
):
    model = CarrierMobilityModel(
        effective_mass=0.26 * SI.me0,
        lattice_mobility_reference=1_500.0,
        impurity_mobility_reference=1_500.0,
        impurity_concentration_reference=1.0e16,
    )

    with pytest.raises(ValueError, match="temperature"):
        model.lattice_mobility(0.0)
