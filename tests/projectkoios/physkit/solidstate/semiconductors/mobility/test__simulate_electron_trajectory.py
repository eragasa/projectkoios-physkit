"""Determinism checks for the pedagogical carrier trajectory."""

import numpy as np

from projectkoios.physkit.constants import SI
from projectkoios.physkit.solidstate.semiconductors.mobility import (
    CarrierMobilityModel,
    simulate_electron_trajectory,
)


def test__simulate_electron_trajectory__fixed_seed_is_reproducible() -> None:
    model = CarrierMobilityModel(
        effective_mass=0.26 * SI.me0,
        lattice_mobility_reference=1_500.0,
        impurity_mobility_reference=1_500.0,
        impurity_concentration_reference=1.0e16,
    )

    first = simulate_electron_trajectory(
        model=model,
        step_count=100,
        seed=42,
    )
    second = simulate_electron_trajectory(
        model=model,
        step_count=100,
        seed=42,
    )

    np.testing.assert_array_equal(first.time, second.time)
    np.testing.assert_array_equal(first.position, second.position)
    np.testing.assert_array_equal(first.velocity, second.velocity)
    assert first.time.shape == (101,)
    assert np.all(np.isfinite(first.position))
