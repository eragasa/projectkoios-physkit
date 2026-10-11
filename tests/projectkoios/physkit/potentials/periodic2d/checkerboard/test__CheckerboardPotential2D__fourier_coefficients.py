"""Fourier tests for ``CheckerboardPotential2D``."""

import numpy as np

from projectkoios.physkit.potentials.periodic2d.checkerboard import (
    CheckerboardPotential2D,
)
from projectkoios.physkit.units import ScalarQuantity, Unitless


def test__fourier_coefficients__couples_four_checkerboard_modes() -> None:
    """The fundamental product has four equal Fourier coefficients."""
    potential = CheckerboardPotential2D(
        amplitude=ScalarQuantity(8.0, Unitless()),
    )
    modes = ((-1, -1), (-1, 1), (0, 0), (1, -1), (1, 1))

    coefficients = potential.fourier_coefficients(
        mode_indices=modes
    ).coefficients.magnitude

    np.testing.assert_allclose(coefficients, np.asarray([2, 2, 0, 2, 2]))
    assert potential.fourier_support == modes[:2] + modes[3:]
