"""Plotting tests for ``BandStructureVisualizer``."""

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from projectkoios.physkit.solidstate.visualization.bands import (
    BandStructurePlotData,
    BandStructureVisualizer,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
    VectorQuantity,
)


def _data() -> BandStructurePlotData:
    return BandStructurePlotData(
        path_coordinates=VectorQuantity(
            magnitude=np.asarray([0.0, 0.5, 1.0]),
            unit=Unitless(),
        ),
        band_energies=MatrixQuantity(
            magnitude=np.asarray(
                [[0.0, 2.0], [0.5, 1.5], [1.0, 1.0]],
            ),
            unit=PhysicalUnit("electron_volt"),
        ),
    )


def test__plot__uses_dimension_independent_path_and_optional_reference() -> None:
    """The visualizer only requires scalar path coordinates and band columns."""
    data = _data()
    reference = MatrixQuantity(
        magnitude=np.asarray(
            [[0.1, 2.1], [0.6, 1.6], [1.1, 1.1]],
        ),
        unit=data.band_energies.unit,
    )
    figure, axes = plt.subplots()

    returned = BandStructureVisualizer(data=data).plot(
        ax=axes,
        reference_band_energies=reference,
        path_label="Normalized path",
        title="Compared bands",
    )

    assert returned is axes
    assert len(axes.lines) == 4
    np.testing.assert_allclose(axes.lines[0].get_xdata(), [0.0, 0.5, 1.0])
    assert axes.get_xlabel() == "Normalized path"
    assert axes.get_ylabel() == "Energy (electron_volt)"
    assert axes.get_title() == "Compared bands"
    plt.close(figure)


def test__plot__rejects_reference_with_mismatched_shape() -> None:
    """Reference columns and path rows must match the primary band matrix."""
    data = _data()
    reference = MatrixQuantity(
        magnitude=np.ones((3, 1)),
        unit=data.band_energies.unit,
    )

    with pytest.raises(ValueError, match="match the primary band shape"):
        BandStructureVisualizer(data=data).plot(
            reference_band_energies=reference
        )
