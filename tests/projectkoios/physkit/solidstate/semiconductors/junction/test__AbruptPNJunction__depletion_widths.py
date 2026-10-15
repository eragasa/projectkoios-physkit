"""Numerical checks for abrupt-junction depletion widths."""

import numpy as np

from projectkoios.physkit.solidstate.semiconductors.junction import (
    AbruptPNJunction,
)


def test__AbruptPNJunction__depletion_widths__preserve_charge_balance(
    silicon_junction: AbruptPNJunction,
) -> None:
    widths = silicon_junction.depletion_widths()

    np.testing.assert_allclose(widths.total, widths.p_side + widths.n_side)
    np.testing.assert_allclose(
        silicon_junction.doping.acceptor_density * widths.p_side,
        silicon_junction.doping.donor_density * widths.n_side,
    )
    np.testing.assert_allclose(
        silicon_junction.built_in_potential(),
        0.794736,
        rtol=1.0e-5,
    )
