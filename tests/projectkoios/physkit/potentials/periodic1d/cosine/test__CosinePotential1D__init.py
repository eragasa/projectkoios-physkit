"""Construction tests for ``CosinePotential1D``."""

import pytest

from projectkoios.physkit.potentials.periodic1d.cosine import (
    CosinePotential1D,
)
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity, Unitless


def test__init__rejects_nonenergy_physical_amplitude() -> None:
    """Physical amplitudes must have energy dimensionality."""
    with pytest.raises(ValueError, match="compatible with energy"):
        CosinePotential1D(
            amplitude=ScalarQuantity(1.0, PhysicalUnit("meter")),
        )


@pytest.mark.parametrize("harmonic", [0, -1])
def test__init__rejects_nonpositive_harmonic(harmonic: int) -> None:
    """Cosine harmonics identify positive reciprocal-mode magnitudes."""
    with pytest.raises(ValueError, match="harmonics must be positive"):
        CosinePotential1D(
            amplitude=ScalarQuantity(1.0, Unitless()),
            harmonic=harmonic,
        )
