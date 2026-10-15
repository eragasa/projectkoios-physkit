"""Verification of ``UnitSystem.time_unit``."""

from projectkoios.physkit.units import PhysicalUnit, Unitless, UnitSystem


def test_returns_selected_time_scale() -> None:
    assert isinstance(UnitSystem.NONDIMENSIONAL.time_unit, Unitless)
    assert UnitSystem.SI.time_unit == PhysicalUnit("second")
    assert UnitSystem.METAL.time_unit == PhysicalUnit("picosecond")
