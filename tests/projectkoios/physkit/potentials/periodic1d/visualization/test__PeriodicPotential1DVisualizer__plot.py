"""Plotting tests for ``PeriodicPotential1DVisualizer``."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from projectkoios.physkit.potentials.periodic1d.cosine import (
    CosinePotential1D,
)
from projectkoios.physkit.potentials.periodic1d.visualization import (
    PeriodicPotential1DVisualizer,
)
from projectkoios.physkit.units import ScalarQuantity, Unitless, VectorQuantity


def test__plot__evaluates_model_on_fractional_coordinates() -> None:
    """The plotted values come from the maintained potential evaluator."""
    potential = CosinePotential1D(
        amplitude=ScalarQuantity(2.0, Unitless()),
    )
    coordinates = VectorQuantity(
        magnitude=np.asarray([0.0, 0.25, 0.5, 1.0]),
        unit=Unitless(),
    )
    figure, axes = plt.subplots()

    returned = PeriodicPotential1DVisualizer(model=potential).plot(
        fractional_coordinates=coordinates,
        ax=axes,
    )

    assert returned is axes
    np.testing.assert_allclose(axes.lines[0].get_xdata(), coordinates.magnitude)
    np.testing.assert_allclose(
        axes.lines[0].get_ydata(),
        np.asarray([2.0, 0.0, -2.0, 2.0]),
        atol=1.0e-14,
    )
    assert axes.get_ylabel() == "Potential (dimensionless)"
    plt.close(figure)
