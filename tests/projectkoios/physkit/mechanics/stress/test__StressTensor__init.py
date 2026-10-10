"""Verification of immutable unit-aware stress-tensor construction."""

import numpy as np
import pytest

from projectkoios.physkit.mechanics.stress import (
    StressSignConvention,
    StressTensor,
)
from projectkoios.physkit.units import MatrixQuantity, PhysicalUnit, Unitless


def test__StressTensor__init__retains_an_immutable_pressure_quantity() -> None:
    source = np.asarray(
        (
            (1.0, 0.1, 0.2),
            (0.3, 2.0, 0.4),
            (0.5, 0.6, 3.0),
        ),
        dtype=np.float64,
    )
    tensor = StressTensor(
        components=MatrixQuantity(source, PhysicalUnit("gigapascal")),
        sign_convention=StressSignConvention.TENSION_POSITIVE,
    )

    source[0, 0] = 99.0

    assert tensor.components.magnitude.tolist() == [
        [1.0, 0.1, 0.2],
        [0.3, 2.0, 0.4],
        [0.5, 0.6, 3.0],
    ]
    assert tensor.components.magnitude.flags.writeable is False
    with pytest.raises(ValueError, match="read-only"):
        tensor.components.magnitude[0, 0] = 4.0


def test__StressTensor__init__rejects_non_three_dimensional_components() -> None:
    with pytest.raises(ValueError, match=r"shape \(3, 3\)"):
        StressTensor(
            components=MatrixQuantity(
                np.eye(2, dtype=np.float64),
                PhysicalUnit("pascal"),
            ),
            sign_convention=StressSignConvention.TENSION_POSITIVE,
        )


def test__StressTensor__init__rejects_non_pressure_units() -> None:
    with pytest.raises(ValueError, match="pressure dimensions"):
        StressTensor(
            components=MatrixQuantity(
                np.eye(3, dtype=np.float64),
                PhysicalUnit("joule"),
            ),
            sign_convention=StressSignConvention.TENSION_POSITIVE,
        )

    with pytest.raises(TypeError, match="physical unit"):
        StressTensor(
            components=MatrixQuantity(
                np.eye(3, dtype=np.float64),
                Unitless(),
            ),
            sign_convention=StressSignConvention.TENSION_POSITIVE,
        )
