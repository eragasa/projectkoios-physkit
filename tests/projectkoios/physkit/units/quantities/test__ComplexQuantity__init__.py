"""Verification of finite complex scalar quantities."""

from dataclasses import FrozenInstanceError

import pytest

from projectkoios.physkit.units.quantities import (
    ComplexQuantity,
    PhysicalUnit,
    Unitless,
)


def test__init__finite_builtin_complex__retains_value_and_unit() -> None:
    """Exact finite complex values retain their nominal unit."""
    quantity = ComplexQuantity(1.0 - 2.0j, PhysicalUnit("electron_volt"))

    assert type(quantity.magnitude) is complex
    assert quantity.magnitude == 1.0 - 2.0j
    assert quantity.unit == PhysicalUnit("electron_volt")


def test__init__noncomplex_or_nonfinite_value__rejects() -> None:
    """Real coercions, booleans, strings, and nonfinite components are rejected."""
    for value in (1.0, 1, True, "1+2j"):
        with pytest.raises(TypeError, match="built-in complex"):
            ComplexQuantity(value, Unitless())  # type: ignore[arg-type]
    for value in (
        complex(float("nan"), 0.0),
        complex(0.0, float("inf")),
    ):
        with pytest.raises(ValueError, match="finite"):
            ComplexQuantity(value, Unitless())


def test__init__invalid_unit__rejects() -> None:
    """Only nominal PhysKit model-system units are admitted."""
    with pytest.raises(TypeError, match="PhysicalUnit or Unitless"):
        ComplexQuantity(1.0 + 0.0j, "electron_volt")  # type: ignore[arg-type]


def test__mutation__frozen_quantity__rejects() -> None:
    """Complex scalar quantities are operationally immutable."""
    quantity = ComplexQuantity(1.0 + 2.0j, Unitless())

    with pytest.raises(FrozenInstanceError):
        quantity.magnitude = 0.0j  # type: ignore[misc]
