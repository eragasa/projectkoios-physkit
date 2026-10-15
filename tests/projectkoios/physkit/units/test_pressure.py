# tests/projectkoios/physkit/units/test_pressure.py

import numpy as np
import pytest

from projectkoios.physkit.units import Pressure
from projectkoios.physkit.units.protocols import UnitQuantityProtocol

# ---------------------------------------------------------------------
# API / protocol compliance
# ---------------------------------------------------------------------


@pytest.mark.unit
def test_pressure_conforms_to_quantity_protocol():
    assert isinstance(Pressure, UnitQuantityProtocol)


@pytest.mark.unit
def test_units_is_intenum():
    from enum import IntEnum

    assert issubclass(Pressure.Units, IntEnum)


# ---------------------------------------------------------------------
# Canonical unit sanity
# ---------------------------------------------------------------------


@pytest.mark.unit
def test_canonical_unit_is_pascal():
    val = Pressure.convert(from_=[1.0, Pressure.Units.Pa], to=Pressure.Units.Pa)
    assert val == 1.0


# ---------------------------------------------------------------------
# Known reference conversions (scalar)
# ---------------------------------------------------------------------


@pytest.mark.unit
def test_atm_to_pa():
    val = Pressure.convert(from_=[1.0, Pressure.Units.atm], to=Pressure.Units.Pa)
    assert abs(val - 101_325.0) < 1e-12


@pytest.mark.unit
def test_pa_to_atm():
    val = Pressure.convert(from_=[101_325.0, Pressure.Units.Pa], to=Pressure.Units.atm)
    assert abs(val - 1.0) < 1e-12


@pytest.mark.unit
def test_psi_to_pa():
    val = Pressure.convert(from_=[1.0, Pressure.Units.psi], to=Pressure.Units.Pa)
    assert abs(val - 6894.757293168) < 1e-9


@pytest.mark.unit
def test_bar_to_pa():
    val = Pressure.convert(from_=[1.0, Pressure.Units.bar], to=Pressure.Units.Pa)
    assert abs(val - 1.0e5) < 1e-12


@pytest.mark.unit
def test_barye_to_pa():
    val = Pressure.convert(from_=[10.0, Pressure.Units.Ba], to=Pressure.Units.Pa)
    assert abs(val - 1.0) < 1e-12


# ---------------------------------------------------------------------
# Round-trip invariants (physics sanity)
# ---------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("unit", list(Pressure.Units))
def test_round_trip_scalar(unit):
    x = 123.456
    y = Pressure.convert(from_=[x, unit], to=Pressure.Units.Pa)
    z = Pressure.convert(from_=[y, Pressure.Units.Pa], to=unit)
    assert abs(z - x) < 1e-12


# ---------------------------------------------------------------------
# Vectorization behavior
# ---------------------------------------------------------------------


@pytest.mark.unit
def test_vectorized_conversion():
    arr = np.array([0.5, 1.0, 2.0])
    out = Pressure.convert(from_=[arr, Pressure.Units.atm], to=Pressure.Units.Pa)
    assert isinstance(out, np.ndarray)
    assert np.allclose(out, arr * 101_325.0)


# ---------------------------------------------------------------------
# Consistency between equivalent units
# ---------------------------------------------------------------------


@pytest.mark.unit
def test_torr_and_mmhg_equivalence():
    x = 760.0
    pa_torr = Pressure.convert(from_=[x, Pressure.Units.Torr], to=Pressure.Units.Pa)
    pa_mmhg = Pressure.convert(from_=[x, Pressure.Units.mmHg], to=Pressure.Units.Pa)
    assert abs(pa_torr - pa_mmhg) < 1e-12


# ---------------------------------------------------------------------
# Error behavior (minimal but explicit)
# ---------------------------------------------------------------------


@pytest.mark.unit
def test_invalid_unit_raises():
    class FakeUnit:
        pass

    with pytest.raises(TypeError):
        Pressure.convert(from_=[1.0, FakeUnit()], to=Pressure.Units.Pa)
