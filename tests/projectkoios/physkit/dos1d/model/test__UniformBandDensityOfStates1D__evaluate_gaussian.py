"""Gaussian DOS tests for ``UniformBandDensityOfStates1D``."""

import numpy as np

from projectkoios.physkit.dos1d.model import (
    UniformBandDensityOfStates1D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    VectorQuantity,
)


def test__evaluate_gaussian__normalizes_one_band_per_cell_length() -> None:
    """A complete spinless 1D band integrates to one state per cell length."""
    model = UniformBandDensityOfStates1D(
        brillouin_zone_measure=ScalarQuantity(
            2.0 * np.pi / 5.0,
            PhysicalUnit("1 / angstrom"),
        ),
    )
    energies = MatrixQuantity(
        magnitude=np.asarray([[0.0], [0.5], [1.0], [1.5]]),
        unit=PhysicalUnit("electron_volt"),
    )

    result = model.evaluate_gaussian(
        band_energies=energies,
        energy_grid=VectorQuantity(
            magnitude=np.linspace(-1.0, 2.5, 3001),
            unit=energies.unit,
        ),
        broadening=ScalarQuantity(0.05, energies.unit),
    )

    assert np.isclose(result.represented_state_density.magnitude, 0.2)
    assert np.isclose(
        result.integrated_state_density.magnitude,
        result.represented_state_density.magnitude,
        rtol=1.0e-10,
    )
    assert result.density_of_states.unit == PhysicalUnit(
        "(1 / angstrom) / (electron_volt)"
    )
