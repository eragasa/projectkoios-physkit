"""Verification of scalar Pint conversion."""

import pytest

from projectkoios.physkit.units.quantities import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
)


class TestPintUnitConverterConvertScalar:
    """Verify multiplicative, affine, and incompatible scalar conversion."""

    def test__convert_scalar_converts_compatible_length(self) -> None:
        result = ModelSystemUnitConverter.convert_scalar(
            quantity=ScalarQuantity(
                magnitude=100.0,
                unit=PhysicalUnit(expression="centimeter"),
            ),
            target=PhysicalUnit(expression="meter"),
        )

        assert result == ScalarQuantity(
            magnitude=1.0,
            unit=PhysicalUnit(expression="meter"),
        )

    def test__convert_scalar_converts_affine_temperature(self) -> None:
        result = ModelSystemUnitConverter.convert_scalar(
            quantity=ScalarQuantity(
                magnitude=26.85,
                unit=PhysicalUnit(expression="degree_Celsius"),
            ),
            target=PhysicalUnit(expression="kelvin"),
        )

        assert result.magnitude == pytest.approx(300.0)
        assert result.unit == PhysicalUnit(expression="kelvin")

    def test__convert_scalar_rejects_incompatible_dimensions(self) -> None:
        with pytest.raises(ValueError, match="dimensionally incompatible"):
            ModelSystemUnitConverter.convert_scalar(
                quantity=ScalarQuantity(
                    magnitude=1.0,
                    unit=PhysicalUnit(expression="meter"),
                ),
                target=PhysicalUnit(expression="kilogram"),
            )
