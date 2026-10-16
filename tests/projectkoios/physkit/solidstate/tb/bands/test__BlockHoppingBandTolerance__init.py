"""Caller-supplied block-band tolerance invariant tests."""

from typing import Any, cast

import pytest

from projectkoios.physkit.solidstate.tb.bands import BlockHoppingBandTolerance
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity


def test__init__accepts_finite_nonnegative_numerical_tolerances() -> None:
    absolute = ScalarQuantity(1.0e-9, PhysicalUnit("electron_volt"))

    tolerance = BlockHoppingBandTolerance(absolute, 1.0e-7)

    assert tolerance.hermiticity_absolute is absolute
    assert tolerance.hermiticity_relative == 1.0e-7


@pytest.mark.parametrize("magnitude", [-1.0, -0.1])
def test__init__rejects_negative_absolute_tolerance(magnitude: float) -> None:
    with pytest.raises(ValueError, match="absolute must be nonnegative"):
        BlockHoppingBandTolerance(
            ScalarQuantity(magnitude, PhysicalUnit("electron_volt")),
            0.0,
        )


def test__init__rejects_wrong_absolute_type() -> None:
    with pytest.raises(TypeError, match="absolute must be ScalarQuantity"):
        BlockHoppingBandTolerance(cast(Any, 0.0), 0.0)


@pytest.mark.parametrize("relative", [-0.1, 1.0, 2.0])
def test__init__rejects_relative_tolerance_outside_half_open_unit_interval(
    relative: float,
) -> None:
    with pytest.raises(ValueError, match=r"lie in \[0, 1\)"):
        BlockHoppingBandTolerance(
            ScalarQuantity(0.0, PhysicalUnit("electron_volt")),
            relative,
        )


def test__init__rejects_wrong_or_nonfinite_relative_tolerance() -> None:
    absolute = ScalarQuantity(0.0, PhysicalUnit("electron_volt"))
    with pytest.raises(TypeError, match="built-in float"):
        BlockHoppingBandTolerance(absolute, cast(Any, 0))
    with pytest.raises(ValueError, match="must be finite"):
        BlockHoppingBandTolerance(absolute, float("nan"))
