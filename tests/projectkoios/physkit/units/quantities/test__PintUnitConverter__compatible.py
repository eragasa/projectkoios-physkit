"""Verification of Pint dimensional compatibility."""

from projectkoios.physkit.units.quantities import (
    ModelSystemUnitConverter,
    PhysicalUnit,
)


def test_reports_dimensional_compatibility() -> None:
    assert ModelSystemUnitConverter.compatible(
        PhysicalUnit("joule * second"),
        PhysicalUnit("kilogram * meter ** 2 / second"),
    )
    assert not ModelSystemUnitConverter.compatible(
        PhysicalUnit("meter"),
        PhysicalUnit("kilogram"),
    )
