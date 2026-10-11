"""Construction tests for ``UniformBandDensityOfStates2D``."""

import pytest

from projectkoios.physkit.dos2d.model import (
    UniformBandDensityOfStates2D,
)
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity


def test__init__rejects_one_dimensional_zone_measure() -> None:
    """A two-dimensional DOS requires reciprocal area rather than length."""
    with pytest.raises(ValueError, match="matching the model dimension"):
        UniformBandDensityOfStates2D(
            brillouin_zone_measure=ScalarQuantity(
                1.0,
                PhysicalUnit("1 / angstrom"),
            )
        )
