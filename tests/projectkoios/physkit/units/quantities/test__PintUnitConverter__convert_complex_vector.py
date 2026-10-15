"""Verification of complex-vector Pint conversion."""

import numpy as np

from projectkoios.physkit.units.quantities import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    ComplexVectorQuantity,
    PhysicalUnit,
)


def test_converts_compatible_complex_vector() -> None:
    converted = MODEL_SYSTEM_UNIT_CONVERTER.convert_complex_vector(
        ComplexVectorQuantity(
            np.array([100.0 + 50.0j]),
            PhysicalUnit("centimeter"),
        ),
        PhysicalUnit("meter"),
    )

    np.testing.assert_allclose(converted.magnitude, [1.0 + 0.5j])
    assert converted.unit == PhysicalUnit("meter")
