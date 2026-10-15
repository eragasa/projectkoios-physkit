# tests/projectkoios/physkit/units/test_force.py

import numpy as np
import pytest

from projectkoios.physkit.units import Force
from projectkoios.physkit.units.protocols import UnitQuantityProtocol


# --- API / protocol compliance
@pytest.mark.unit
def test_force_conforms_to_quantity_protocol():
    assert isinstance(Force, UnitQuantityProtocol)


@pytest.mark.unit
def test_units_is_intenum():
    from enum import IntEnum

    assert issubclass(Force.Units, IntEnum)


# Canonical unit sanity
@pytest.mark.unit
def test_canonical_unit_is_newton():
    val = Force.convert(from_=[1.0, Force.Units.N], to=Force.Units.N)
    assert val == 1.0


# --- Known reference conversions (scalar)
@pytest.mark.unit
def test_lbf_to_newton():
    val = Force.convert(from_=[1.0, Force.Units.lbf], to=Force.Units.N)
    assert abs(val - 4.4482216152605) < 1e-12


@pytest.mark.unit
def test_newton_to_lbf():
    val = Force.convert(from_=[4.4482216152605, Force.Units.N], to=Force.Units.lbf)
    assert abs(val - 1.0) < 1e-12


@pytest.mark.unit
def test_dyne_to_newton():
    val = Force.convert(from_=[1.0, Force.Units.dyn], to=Force.Units.N)
    assert abs(val - 1.0e-5) < 1e-15


@pytest.mark.unit
def test_newton_to_dyne():
    val = Force.convert(from_=[1.0, Force.Units.N], to=Force.Units.dyn)
    assert abs(val - 1.0e5) < 1e-9


# Round-trip invariants


@pytest.mark.unit
@pytest.mark.parametrize("unit", list(Force.Units))
def test_round_trip_scalar(unit):
    x = 123.456
    y = Force.convert(from_=[x, unit], to=Force.Units.N)
    z = Force.convert(from_=[y, Force.Units.N], to=unit)
    assert abs(z - x) < 1e-12


# ---Vectorization behavior
@pytest.mark.unit
def test_vectorized_conversion():
    arr = np.array([0.5, 1.0, 2.0])
    out = Force.convert(from_=[arr, Force.Units.lbf], to=Force.Units.N)
    assert isinstance(out, np.ndarray)
    assert np.allclose(out, arr * 4.4482216152605)


# Error behavior


@pytest.mark.unit
def test_invalid_unit_raises():
    class FakeUnit:
        pass

    with pytest.raises(TypeError):
        Force.convert(from_=[1.0, FakeUnit()], to=Force.Units.N)
