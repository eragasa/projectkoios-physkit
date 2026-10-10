"""Verification of explicit stress-unit and sign conversion."""

import numpy as np

from projectkoios.physkit.mechanics.stress import (
    StressSignConvention,
    StressTensor,
)
from projectkoios.physkit.units import MatrixQuantity, PhysicalUnit


def test__StressTensor__to_si__converts_components_to_pascals() -> None:
    tensor = StressTensor(
        components=MatrixQuantity(
            np.diag(np.asarray((1.0, 2.0, 3.0), dtype=np.float64)),
            PhysicalUnit("gigapascal"),
        ),
        sign_convention=StressSignConvention.TENSION_POSITIVE,
    )

    converted = tensor.to_si()

    assert converted.components.unit == PhysicalUnit("pascal")
    assert converted.components.magnitude.tolist() == [
        [1.0e9, 0.0, 0.0],
        [0.0, 2.0e9, 0.0],
        [0.0, 0.0, 3.0e9],
    ]
    assert converted.sign_convention is StressSignConvention.TENSION_POSITIVE


def test__StressTensor__to_sign_convention__reverses_every_component() -> None:
    tensor = StressTensor(
        components=MatrixQuantity(
            np.asarray(
                ((1.0, 0.1, 0.2), (0.3, 2.0, 0.4), (0.5, 0.6, 3.0)),
                dtype=np.float64,
            ),
            PhysicalUnit("kilobar"),
        ),
        sign_convention=StressSignConvention.COMPRESSION_POSITIVE,
    )

    converted = tensor.to_sign_convention(
        StressSignConvention.TENSION_POSITIVE,
    )

    assert converted.components.unit == PhysicalUnit("kilobar")
    assert np.array_equal(
        converted.components.magnitude,
        -tensor.components.magnitude,
    )


def test__StressTensor__hydrostatic_pressure__is_positive_under_compression() -> None:
    compression_positive = StressTensor(
        components=MatrixQuantity(
            np.diag(np.asarray((1.0, 2.0, 3.0), dtype=np.float64)),
            PhysicalUnit("gigapascal"),
        ),
        sign_convention=StressSignConvention.COMPRESSION_POSITIVE,
    )
    tension_positive = compression_positive.to_sign_convention(
        StressSignConvention.TENSION_POSITIVE,
    )

    assert compression_positive.hydrostatic_pressure.magnitude == 2.0
    assert tension_positive.hydrostatic_pressure.magnitude == 2.0
    assert compression_positive.hydrostatic_pressure.unit == PhysicalUnit("gigapascal")
