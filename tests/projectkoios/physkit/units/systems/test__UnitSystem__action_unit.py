"""Verification of ``UnitSystem.action_unit``."""

from projectkoios.physkit.units import PhysicalUnit, Unitless, UnitSystem


def test_returns_selected_action_scale() -> None:
    assert isinstance(UnitSystem.NONDIMENSIONAL.action_unit, Unitless)
    assert UnitSystem.SI.action_unit == PhysicalUnit("joule * second")
    assert UnitSystem.METAL.action_unit == PhysicalUnit("electron_volt * picosecond")
