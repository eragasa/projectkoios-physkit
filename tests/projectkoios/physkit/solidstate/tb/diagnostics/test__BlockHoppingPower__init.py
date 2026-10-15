"""Invariant tests for one unit-aware block power."""

from typing import Any, cast

import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.solidstate.tb.diagnostics import BlockHoppingPower
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity


def test__init__retains_displacement_and_frobenius_power() -> None:
    displacement = LatticeDisplacement(LatticeDimension.ONE, (0,))
    power = ScalarQuantity(3.0, PhysicalUnit("electron_volt ** 2"))

    result = BlockHoppingPower(displacement, power)

    assert result.displacement is displacement
    assert result.frobenius_power is power


def test__init__rejects_wrong_types_and_negative_power() -> None:
    displacement = LatticeDisplacement(LatticeDimension.ONE, (0,))

    with pytest.raises(TypeError, match="displacement must be LatticeDisplacement"):
        BlockHoppingPower(
            cast(Any, object()),
            ScalarQuantity(1.0, PhysicalUnit("electron_volt ** 2")),
        )
    with pytest.raises(TypeError, match="frobenius_power must be ScalarQuantity"):
        BlockHoppingPower(displacement, cast(Any, 1.0))
    with pytest.raises(ValueError, match="must be nonnegative"):
        BlockHoppingPower(
            displacement,
            ScalarQuantity(-1.0, PhysicalUnit("electron_volt ** 2")),
        )
