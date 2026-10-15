"""Verification of ``UnitSystem.length_unit``."""

from projectkoios.physkit.units import PhysicalUnit, Unitless, UnitSystem


def test_returns_selected_length_scale() -> None:
    assert isinstance(UnitSystem.NONDIMENSIONAL.length_unit, Unitless)
    assert UnitSystem.SI.length_unit == PhysicalUnit("meter")
    assert UnitSystem.METAL.length_unit == PhysicalUnit("angstrom")
