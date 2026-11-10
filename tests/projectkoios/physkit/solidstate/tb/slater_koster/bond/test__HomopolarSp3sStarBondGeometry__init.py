"""Tests for ordinary-cell Slater--Koster bond geometry."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.solidstate.tb.slater_koster.bond import (
    HomopolarSp3sStarBondGeometry,
    SlaterKosterDirectedBond,
)
from projectkoios.physkit.units import PhysicalUnit, UnitSystem, VectorQuantity


class TestHomopolarSp3sStarBondGeometry:
    """Keep cell displacement separate from physical center-to-center geometry."""

    def test__init__derives_direction_without_changing_directed_bond(self) -> None:
        directed_bond = SlaterKosterDirectedBond(
            ordinary_cell_displacement=LatticeDisplacement(
                dimension=LatticeDimension.THREE,
                components=(1, -2, 3),
            ),
            source_center_identifier="origin_center",
            target_center_identifier="quarter_diagonal_center",
        )
        geometry = HomopolarSp3sStarBondGeometry(
            directed_bond=directed_bond,
            cartesian_displacement=VectorQuantity(
                magnitude=np.asarray((2.0, -1.0, 2.0)),
                unit=PhysicalUnit(expression="angstrom"),
            ),
            unit_system=UnitSystem.METAL,
        )

        np.testing.assert_array_equal(
            geometry.bond_direction.components.magnitude,
            np.asarray((2.0, -1.0, 2.0)) / 3.0,
        )
        assert geometry.directed_bond is directed_bond
        assert directed_bond.reversed() == SlaterKosterDirectedBond(
            ordinary_cell_displacement=LatticeDisplacement(
                dimension=LatticeDimension.THREE,
                components=(-1, 2, -3),
            ),
            source_center_identifier="quarter_diagonal_center",
            target_center_identifier="origin_center",
        )

    def test__init__uses_atomic_length_by_default(self) -> None:
        geometry = HomopolarSp3sStarBondGeometry(
            directed_bond=SlaterKosterDirectedBond(
                ordinary_cell_displacement=LatticeDisplacement(
                    dimension=LatticeDimension.THREE,
                    components=(0, 0, 0),
                ),
                source_center_identifier="origin_center",
                target_center_identifier="quarter_diagonal_center",
            ),
            cartesian_displacement=VectorQuantity(
                magnitude=np.asarray((1.0, 1.0, 1.0)),
                unit=UnitSystem.ATOMIC.length_unit,
            ),
        )

        assert geometry.unit_system is UnitSystem.ATOMIC
        assert geometry.cartesian_displacement.unit == UnitSystem.ATOMIC.length_unit

    def test__init__rejects_zero_cell_self_block(self) -> None:
        with pytest.raises(ValueError, match="onsite ownership"):
            SlaterKosterDirectedBond(
                ordinary_cell_displacement=LatticeDisplacement(
                    dimension=LatticeDimension.THREE,
                    components=(0, 0, 0),
                ),
                source_center_identifier="origin_center",
                target_center_identifier="origin_center",
            )
