"""Verify coherent atomic-unit ownership."""

import pytest

from projectkoios.physkit.units import PhysicalUnit, UnitSystem

pytestmark = pytest.mark.unit


class TestUnitSystem:
    """Check the exact atomic units and fundamental electron quantities."""

    def test__UnitSystem__atomic_units_are_coherent(self) -> None:
        """Expose one coherent atomic length, mass, time, energy, and action set."""
        assert (
            UnitSystem.ATOMIC.length_unit,
            UnitSystem.ATOMIC.mass_unit,
            UnitSystem.ATOMIC.time_unit,
            UnitSystem.ATOMIC.energy_unit,
            UnitSystem.ATOMIC.action_unit,
        ) == (
            PhysicalUnit("bohr"),
            PhysicalUnit("electron_mass"),
            PhysicalUnit("atomic_unit_of_time"),
            PhysicalUnit("hartree"),
            PhysicalUnit("atomic_unit_of_action"),
        )

    def test__UnitSystem__atomic_hbar_is_unity(self) -> None:
        """Represent reduced Planck's constant as one atomic action unit."""
        assert UnitSystem.ATOMIC.hbar.magnitude == pytest.approx(
            1.0,
            rel=0.0,
            abs=2.0e-15,
        )
