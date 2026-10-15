"""Shared abrupt-junction test constructions."""

import pytest

from projectkoios.physkit.solidstate.semiconductors.junction import (
    AbruptDopingProfile,
    AbruptPNJunction,
    SemiconductorMaterial,
)


@pytest.fixture
def silicon_junction() -> AbruptPNJunction:
    """Return the silicon-like configuration extracted from the notebook."""
    return AbruptPNJunction(
        material=SemiconductorMaterial(
            band_gap=1.12,
            electron_affinity=4.05,
            relative_permittivity=11.7,
            conduction_density_of_states=2.8e25,
            valence_density_of_states=1.04e25,
            electron_mobility=1_350.0e-4,
            hole_mobility=480.0e-4,
            temperature=300.0,
        ),
        doping=AbruptDopingProfile(
            acceptor_density=1.0e23,
            donor_density=1.0e22,
        ),
    )
