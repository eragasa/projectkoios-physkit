"""Numerical verification for the empirical mobility contributions."""

import numpy as np

from projectkoios.physkit.constants import SI
from projectkoios.physkit.solidstate.semiconductors.mobility import (
    CarrierMobilityModel,
)


def _model() -> CarrierMobilityModel:
    return CarrierMobilityModel(
        effective_mass=0.26 * SI.me0,
        lattice_mobility_reference=1_500.0,
        impurity_mobility_reference=1_500.0,
        impurity_concentration_reference=1.0e16,
    )


def test__CarrierMobilityModel__mobilities__match_source_power_laws() -> None:
    model = _model()

    lattice = model.lattice_mobility(np.array([300.0, 600.0]))
    impurity = model.impurity_mobility(
        np.array([300.0, 600.0]),
        np.array([1.0e16, 2.0e16]),
    )
    total = model.total_mobility(300.0, 1.0e16)

    np.testing.assert_allclose(
        lattice,
        np.array([1_500.0, 1_500.0 * 2.0**-1.5]),
    )
    np.testing.assert_allclose(
        impurity,
        np.array([1_500.0, 1_500.0 * 2.0**1.5 / 2.0]),
    )
    np.testing.assert_allclose(total, 750.0)


def test__CarrierMobilityModel__relaxation_time__round_trips_mobility() -> None:
    model = _model()
    mobility = np.array([450.0, 1_350.0])

    reconstructed = model.mobility_from_relaxation_time(model.relaxation_time(mobility))

    np.testing.assert_allclose(reconstructed, mobility)


def test__CarrierMobilityModel__thermal_velocity__matches_equipartition() -> None:
    model = _model()

    actual = model.thermal_velocity(300.0)
    expected = np.sqrt(3.0 * SI.k_B * 300.0 / model.effective_mass)

    np.testing.assert_allclose(actual, expected)
