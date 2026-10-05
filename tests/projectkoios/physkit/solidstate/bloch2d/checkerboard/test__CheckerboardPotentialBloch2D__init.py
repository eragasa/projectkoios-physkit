"""Checkerboard model construction requirements."""

import numpy as np
import pytest

from projectkoios.physkit.solidstate.bloch2d.checkerboard.model import (
    CheckerboardPotentialBloch2D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.cell import (
    PeriodicPrimitiveCell2D,
)
from projectkoios.physkit.units import MatrixQuantity, ScalarQuantity, UnitSystem


def test__init__rejects_a_nonsquare_cell() -> None:
    cell = PeriodicPrimitiveCell2D(
        primitive_basis=MatrixQuantity(
            magnitude=np.diag([8.0, 10.0]),
            unit=UnitSystem.METAL.length_unit,
        ),
        unit_system=UnitSystem.METAL,
    )

    with pytest.raises(
        ValueError,
        match="cell must be square with orthogonal equal-length vectors",
    ):
        CheckerboardPotentialBloch2D(
            cell=cell,
            amplitude=ScalarQuantity(1.0, UnitSystem.METAL.energy_unit),
        )
