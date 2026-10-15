"""Compatibility tests for ``CosinePotentialFourier1D``."""

import numpy as np

from projectkoios.physkit.periodic.lattice import DirectLattice1D
from projectkoios.physkit.solidstate.electronic.models import (
    CosinePotentialFourier1D,
)


def test__coefficients__delegates_integer_mode_couplings() -> None:
    """The wavevector API retains its two fundamental cosine coefficients."""
    lattice = DirectLattice1D(lattice_constant=2.0)
    potential = CosinePotentialFourier1D(amplitude=6.0, lattice=lattice)
    fundamental = np.pi

    coefficients = potential.coefficients(
        np.asarray(
            [-2.0 * fundamental, -fundamental, 0.5, fundamental, 2.0 * fundamental]
        )
    )

    np.testing.assert_allclose(coefficients, np.asarray([0.0, 3.0, 0.0, 3.0, 0.0]))
