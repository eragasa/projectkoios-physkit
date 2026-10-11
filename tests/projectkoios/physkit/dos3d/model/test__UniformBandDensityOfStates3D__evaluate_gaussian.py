"""Normalization tests for ``UniformBandDensityOfStates3D``."""

import numpy as np

from projectkoios.physkit.dos3d.model import (
    UniformBandDensityOfStates3D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


def test__evaluate_gaussian__normalizes_one_spin_degenerate_band() -> None:
    """One complete band with spin two represents two states per unit volume."""
    model = UniformBandDensityOfStates3D(
        brillouin_zone_measure=ScalarQuantity((2.0 * np.pi) ** 3, Unitless()),
        state_degeneracy=2.0,
    )
    result = model.evaluate_gaussian(
        band_energies=MatrixQuantity(
            magnitude=np.asarray([[0.0], [0.5], [1.0]]),
            unit=Unitless(),
        ),
        energy_grid=VectorQuantity(
            magnitude=np.linspace(-2.0, 3.0, 2001),
            unit=Unitless(),
        ),
        broadening=ScalarQuantity(0.08, Unitless()),
    )

    assert result.represented_state_density.magnitude == 2.0
    assert np.isclose(result.integrated_state_density.magnitude, 2.0)
