"""Tests for explicit Slater--Koster bond directions."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.solidstate.tb.slater_koster.bond import (
    SlaterKosterBondDirection,
)
from projectkoios.physkit.units import PhysicalUnit, UnitSystem, VectorQuantity


class TestSlaterKosterBondDirection:
    """Require an explicit three-component nondimensional unit vector."""

    def test__init__rejects_invalid_direction_quantities(self) -> None:
        with pytest.raises(ValueError, match="unit Euclidean norm"):
            SlaterKosterBondDirection(
                components=VectorQuantity(
                    magnitude=np.asarray((1.0, 1.0, 1.0)),
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                )
            )
        with pytest.raises(ValueError, match="must be nondimensional"):
            SlaterKosterBondDirection(
                components=VectorQuantity(
                    magnitude=np.asarray((1.0, 0.0, 0.0)),
                    unit=PhysicalUnit(expression="angstrom"),
                )
            )
        with pytest.raises(TypeError, match="excluding booleans"):
            VectorQuantity(
                magnitude=np.asarray((True, False, False)),
                unit=UnitSystem.NONDIMENSIONAL.length_unit,
            )
