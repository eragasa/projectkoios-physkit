"""Plotting tests for ``DensityOfStatesVisualizer``."""

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from projectkoios.physkit.dos.visualization import (
    DensityOfStatesPlotData,
    DensityOfStatesVisualizer,
)
from projectkoios.physkit.units import PhysicalUnit, VectorQuantity


def test__plot__supports_vertical_energy_for_band_composites() -> None:
    """Vertical energy orientation is independent of reciprocal-space dimension."""
    data = DensityOfStatesPlotData(
        energies=VectorQuantity(
            magnitude=np.asarray([0.0, 1.0, 2.0]),
            unit=PhysicalUnit("electron_volt"),
        ),
        density_of_states=VectorQuantity(
            magnitude=np.asarray([0.0, 3.0, 1.0]),
            unit=PhysicalUnit("1 / (electron_volt * angstrom ** 2)"),
        ),
    )
    figure, axes = plt.subplots()

    returned = DensityOfStatesVisualizer(data=data).plot(
        ax=axes,
        orientation="energy_vertical",
        label="two-dimensional DOS",
    )

    assert returned is axes
    np.testing.assert_allclose(axes.lines[0].get_xdata(), [0.0, 3.0, 1.0])
    np.testing.assert_allclose(axes.lines[0].get_ydata(), [0.0, 1.0, 2.0])
    assert axes.lines[0].get_label() == "two-dimensional DOS"
    plt.close(figure)
