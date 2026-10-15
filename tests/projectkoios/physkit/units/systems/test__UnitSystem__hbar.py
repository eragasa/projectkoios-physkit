"""Verification of ``UnitSystem.hbar``."""

import numpy as np

from projectkoios.physkit.units import PhysicalUnit, Unitless, UnitSystem


def test_returns_reduced_planck_constant_in_selected_scale() -> None:
    assert UnitSystem.NONDIMENSIONAL.hbar.magnitude == 1.0
    assert isinstance(UnitSystem.NONDIMENSIONAL.hbar.unit, Unitless)
    assert UnitSystem.SI.hbar.unit == PhysicalUnit("joule * second")
    assert UnitSystem.METAL.hbar.unit == PhysicalUnit("electron_volt * picosecond")
    np.testing.assert_allclose(
        UnitSystem.METAL.hbar.magnitude,
        6.582_119_565_476_075e-4,
        rtol=5.0e-10,
        atol=0.0,
    )
