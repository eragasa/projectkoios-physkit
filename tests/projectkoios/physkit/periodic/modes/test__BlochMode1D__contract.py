"""Contract tests for Bloch modes on the canonical reciprocal basis."""

from __future__ import annotations

import numpy as np

from projectkoios.physkit.periodic.basis import ReciprocalModeBasis1D
from projectkoios.physkit.periodic.lattice import ReciprocalLattice1D
from projectkoios.physkit.periodic.modes import BlochMode1D


def test__values__reconstructs_cell_periodic_and_carrier_factors() -> None:
    """Bloch reconstruction consumes vectors from the canonical basis chain."""
    reciprocal = ReciprocalLattice1D(2.0)
    basis = ReciprocalModeBasis1D(reciprocal, [-1, 0, 1])
    coefficients = np.asarray([1.0 + 0.0j, 2.0 - 1.0j, -0.5j])
    mode = BlochMode1D(0.25, basis, coefficients)
    coordinates = np.asarray([-0.2, 0.0, 0.3])
    plane_waves = np.exp(1j * basis.vectors[:, None] * coordinates[None, :])
    expected_periodic = coefficients @ plane_waves

    np.testing.assert_allclose(mode.cell_periodic_part(coordinates), expected_periodic)
    np.testing.assert_allclose(
        mode.values(coordinates),
        np.exp(1j * mode.wavevector * coordinates) * expected_periodic,
    )
    np.testing.assert_allclose(mode.reciprocal_weights(), np.abs(coefficients) ** 2)
