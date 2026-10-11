"""Construction tests for ``BandStructurePlotData``."""

import numpy as np
import pytest

from projectkoios.physkit.solidstate.visualization.bands import (
    BandStructurePlotData,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
    VectorQuantity,
)


def test__init__rejects_nonenergy_band_units() -> None:
    """Band ordinates must represent energy or an explicit nondimensional scale."""
    with pytest.raises(ValueError, match="compatible with energy"):
        BandStructurePlotData(
            path_coordinates=VectorQuantity(
                magnitude=np.asarray([0.0, 1.0]),
                unit=Unitless(),
            ),
            band_energies=MatrixQuantity(
                magnitude=np.ones((2, 1)),
                unit=PhysicalUnit("meter"),
            ),
        )
