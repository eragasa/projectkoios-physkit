"""Numerical checks for the abrupt-junction equilibrium evaluation."""

import numpy as np

from projectkoios.physkit.solidstate.semiconductors.junction import (
    AbruptPNJunction,
)


def test__AbruptPNJunction__solve_equilibrium__returns_correlated_finite_profiles(
    silicon_junction: AbruptPNJunction,
) -> None:
    result = silicon_junction.solve_equilibrium(point_count=101)
    intrinsic_density = silicon_junction.material.intrinsic_carrier_density

    assert result.position.shape == (101,)
    np.testing.assert_allclose(result.electrostatic_potential[0], 0.0, atol=1e-14)
    np.testing.assert_allclose(
        result.electrostatic_potential[-1],
        silicon_junction.built_in_potential(),
    )
    np.testing.assert_allclose(
        result.carriers.electron * result.carriers.hole,
        intrinsic_density**2,
        rtol=1e-14,
    )
    np.testing.assert_allclose(
        result.currents.electron_total,
        result.currents.electron_drift + result.currents.electron_diffusion,
    )
    np.testing.assert_allclose(
        result.currents.hole_total,
        result.currents.hole_drift + result.currents.hole_diffusion,
    )
    assert np.all(np.isfinite(result.electric_field))
    assert np.min(result.charge_density) < 0.0
    assert np.max(result.charge_density) > 0.0
