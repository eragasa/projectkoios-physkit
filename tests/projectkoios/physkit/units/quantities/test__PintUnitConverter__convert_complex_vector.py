"""Verification of complex-vector Pint conversion."""

import numpy as np

from projectkoios.physkit.units.quantities import (
    ComplexVectorQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
)


def test_converts_compatible_complex_vector() -> None:
    converted = ModelSystemUnitConverter.convert_complex_vector(
        ComplexVectorQuantity(
            np.array([100.0 + 50.0j]),
            PhysicalUnit("centimeter"),
        ),
        PhysicalUnit("meter"),
    )

    np.testing.assert_allclose(converted.magnitude, [1.0 + 0.5j])
    assert converted.unit == PhysicalUnit("meter")
