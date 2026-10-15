"""Real-space tests for ``CheckerboardPotential2D``."""

import numpy as np

from projectkoios.physkit.potentials.periodic2d.checkerboard import (
    CheckerboardPotential2D,
)
from projectkoios.physkit.units import MatrixQuantity, ScalarQuantity, Unitless


def test__evaluate__matches_separable_checkerboard_values() -> None:
    """The two-dimensional model is the product of its cosine factors."""
    potential = CheckerboardPotential2D(
        amplitude=ScalarQuantity(8.0, Unitless()),
    )

    values = potential.evaluate(
        fractional_coordinates=MatrixQuantity(
            magnitude=np.asarray([[0.0, 0.0], [0.5, 0.0], [0.25, 0.25], [1.0, 1.0]]),
            unit=Unitless(),
        )
    ).potential_values

    np.testing.assert_allclose(
        values.magnitude,
        np.asarray([8.0, -8.0, 0.0, 8.0]),
        atol=1.0e-14,
    )
