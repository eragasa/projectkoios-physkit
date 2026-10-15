"""Verification of ``UnitSystem.mass_unit``."""

from projectkoios.physkit.units import PhysicalUnit, Unitless, UnitSystem


def test_returns_selected_particle_mass_scale() -> None:
    assert isinstance(UnitSystem.NONDIMENSIONAL.mass_unit, Unitless)
    assert UnitSystem.SI.mass_unit == PhysicalUnit("kilogram")
    assert UnitSystem.METAL.mass_unit == PhysicalUnit("dalton")
