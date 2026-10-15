# tests/projectkoios/physkit/units/test_temperature.py
# Eugene Joseph M. Ragasa, 2026
import inspect
from enum import IntEnum

import numpy as np

from projectkoios.physkit.numeric import as_f64_array
from projectkoios.physkit.units import Temperature
from projectkoios.physkit.units.protocols import UnitQuantityProtocol


# --- UnitQuantityProtocol conformity test ---
def test_Temperature_is_UnitQuantityProtocol():
    assert isinstance(Temperature, UnitQuantityProtocol)


def test_Temperature_has_Units_enum():
    assert hasattr(Temperature, "Units")


def test_Temperature_Units_is_IntEnum():
    assert issubclass(Temperature.Units, IntEnum)


def test_Temperature_convert_exists():
    assert hasattr(Temperature, "convert")
    assert callable(Temperature.convert)


def test_Temperature_convert_check_signature():
    sig_Temperature_convert = inspect.signature(Temperature.convert)
    parameters_Temperature_convert = list(sig_Temperature_convert.parameters.keys())

    sig_UnitQuantityProtocol_convert = inspect.signature(UnitQuantityProtocol.convert)
    parameters_UnitQuantityProtocol_convert = list(
        sig_UnitQuantityProtocol_convert.parameters.keys()
    )

    assert parameters_Temperature_convert == parameters_UnitQuantityProtocol_convert


def test_Temperature_conforms_to_protocol():

    Q = Temperature

    # protocol conformance
    assert isinstance(Q, UnitQuantityProtocol)

    # has Units enum
    assert hasattr(Q, "Units")
    assert issubclass(Q.Units, IntEnum)

    # has convert
    assert callable(Q.convert)

    units = list(Q.Units)
    u0, u1 = units[0], units[1]

    samples = [
        1.0,
        [1.0, 2.0, 3.0],
        [[1.0]],
    ]

    # shape preservation
    for x in samples:
        x_arr = as_f64_array(x)
        y = Q.convert(x, u0, u1)
        y_arr = as_f64_array(y)

        assert y_arr.shape == x_arr.shape

    # round-trip
    for x in samples:
        x_arr = as_f64_array(x)

        y = Q.convert(x, u0, u1)
        x2 = Q.convert(y, u1, u0)
        x2_arr = as_f64_array(x2)

        assert np.allclose(x_arr, x2_arr, atol=1e-12)
