"""Fourier tests for ``CosineProductPotential3D``."""

import numpy as np

from projectkoios.physkit.potentials.periodic3d.cosine_product import (
    CosineProductPotential3D,
)
from projectkoios.physkit.units import ScalarQuantity, Unitless


def test__fourier_coefficients__has_eight_mode_support() -> None:
    """The fundamental three-cosine product has eight equal coefficients."""
    potential = CosineProductPotential3D(
        amplitude=ScalarQuantity(8.0, Unitless()),
    )

    coefficients = potential.fourier_coefficients(
        mode_indices=potential.fourier_support + ((0, 0, 0),)
    ).coefficients.magnitude

    np.testing.assert_allclose(coefficients[:8], np.ones(8))
    assert coefficients[8] == 0.0
