"""Normalization tests for ``UniformBandDensityOfStates2D``."""

import numpy as np

from projectkoios.physkit.dos2d.model import (
    UniformBandDensityOfStates2D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


def test__evaluate_gaussian__normalizes_spin_degenerate_two_band_states() -> None:
    """Two complete bands with spin two represent four states per unit area."""
    model = UniformBandDensityOfStates2D(
        brillouin_zone_measure=ScalarQuantity((2.0 * np.pi) ** 2, Unitless()),
        state_degeneracy=2.0,
    )
    result = model.evaluate_gaussian(
        band_energies=MatrixQuantity(
            magnitude=np.asarray([[0.0, 1.0], [0.5, 1.5]]),
            unit=Unitless(),
        ),
        energy_grid=VectorQuantity(
            magnitude=np.linspace(-2.0, 3.0, 2001),
            unit=Unitless(),
        ),
        broadening=ScalarQuantity(0.08, Unitless()),
    )

    assert result.represented_state_density.magnitude == 4.0
    assert np.isclose(result.integrated_state_density.magnitude, 4.0)
