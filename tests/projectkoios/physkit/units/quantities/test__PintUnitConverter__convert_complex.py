"""Verification of complex-scalar Pint conversion."""

import pytest

from projectkoios.physkit.units.quantities import (
    ComplexQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
)


def test__convert_complex__compatible_units__converts_components() -> None:
    """One multiplicative factor converts both complex components."""
    converted = ModelSystemUnitConverter.convert_complex(
        ComplexQuantity(100.0 + 50.0j, PhysicalUnit("centimeter")),
        PhysicalUnit("meter"),
    )

    assert converted == ComplexQuantity(1.0 + 0.5j, PhysicalUnit("meter"))


def test__convert_complex__affine_or_incompatible_units__rejects() -> None:
    """Complex conversion rejects affine offsets and incompatible dimensions."""
    with pytest.raises(ValueError, match="zero affine offset"):
        ModelSystemUnitConverter.convert_complex(
            ComplexQuantity(1.0 + 2.0j, PhysicalUnit("degree_Celsius")),
            PhysicalUnit("kelvin"),
        )
    with pytest.raises(ValueError, match="dimensionally incompatible"):
        ModelSystemUnitConverter.convert_complex(
            ComplexQuantity(1.0 + 2.0j, PhysicalUnit("meter")),
            PhysicalUnit("second"),
        )


def test__convert_complex__wrong_quantity_type__rejects() -> None:
    """The conversion boundary requires exact nominal complex-quantity ownership."""
    with pytest.raises(TypeError, match="ComplexQuantity"):
        ModelSystemUnitConverter.convert_complex(  # type: ignore[arg-type]
            1.0 + 2.0j,
            PhysicalUnit("meter"),
        )
