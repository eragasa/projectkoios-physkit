"""Evaluation tests for ``CosineProductPotential3D``."""

import numpy as np

from projectkoios.physkit.potentials.periodic3d.cosine_product import (
    CosineProductPotential3D,
)
from projectkoios.physkit.units import MatrixQuantity, ScalarQuantity, Unitless


def test__evaluate__matches_three_dimensional_cosine_product() -> None:
    """The three-dimensional model multiplies all three periodic factors."""
    potential = CosineProductPotential3D(
        amplitude=ScalarQuantity(8.0, Unitless()),
    )

    values = potential.evaluate(
        fractional_coordinates=MatrixQuantity(
            magnitude=np.asarray([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 1.0, 1.0]]),
            unit=Unitless(),
        )
    ).potential_values

    np.testing.assert_allclose(values.magnitude, np.asarray([8.0, -8.0, 8.0]))
