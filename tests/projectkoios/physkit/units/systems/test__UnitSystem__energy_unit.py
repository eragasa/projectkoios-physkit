"""Verification of ``UnitSystem.energy_unit``."""

from projectkoios.physkit.units import PhysicalUnit, Unitless, UnitSystem


def test_returns_selected_energy_scale() -> None:
    assert isinstance(UnitSystem.NONDIMENSIONAL.energy_unit, Unitless)
    assert UnitSystem.SI.energy_unit == PhysicalUnit("joule")
    assert UnitSystem.METAL.energy_unit == PhysicalUnit("electron_volt")
