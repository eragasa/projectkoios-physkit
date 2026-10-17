"""Verification of vector Pint conversion."""

import numpy as np

from projectkoios.physkit.units.quantities import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    VectorQuantity,
)


class TestPintUnitConverterConvertVector:
    """Verify multiplicative and affine vector conversion."""

    def test__convert_vector_converts_affine_temperatures(self) -> None:
        result = ModelSystemUnitConverter.convert_vector(
            quantity=VectorQuantity(
                magnitude=np.array([0.0, 26.85, 100.0], dtype=np.float64),
                unit=PhysicalUnit(expression="degree_Celsius"),
            ),
            target=PhysicalUnit(expression="kelvin"),
        )

        np.testing.assert_allclose(
            result.magnitude,
            np.array([273.15, 300.0, 373.15], dtype=np.float64),
            rtol=0.0,
            atol=6.0e-14,
        )
        assert result.unit == PhysicalUnit(expression="kelvin")
