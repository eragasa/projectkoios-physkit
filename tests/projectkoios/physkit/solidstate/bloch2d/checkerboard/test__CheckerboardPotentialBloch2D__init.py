"""Checkerboard model construction requirements."""

import numpy as np
import pytest

from projectkoios.physkit.solidstate.bloch2d.cell import BlochCell2D
from projectkoios.physkit.solidstate.bloch2d.checkerboard.model import (
    CheckerboardPotentialBloch2D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ScalarQuantity,
    Unitless,
    UnitSystem,
)


def _cell(*, basis: np.ndarray, unit_system: UnitSystem) -> BlochCell2D:
    return BlochCell2D(
        direct_basis=MatrixQuantity(
            magnitude=basis,
            unit=unit_system.length_unit,
        ),
        unit_system=unit_system,
    )


@pytest.mark.parametrize(
    "unit_system",
    [UnitSystem.METAL, UnitSystem.NONDIMENSIONAL],
)
def test__init__accepts_a_rotated_square_in_supported_unit_systems(
    unit_system: UnitSystem,
) -> None:
    angle = np.pi / 7.0
    rotated_square = 3.0 * np.array(
        [
            [np.cos(angle), -np.sin(angle)],
            [np.sin(angle), np.cos(angle)],
        ]
    )

    model = CheckerboardPotentialBloch2D(
        cell=_cell(basis=rotated_square, unit_system=unit_system),
        amplitude=ScalarQuantity(0.8, unit_system.energy_unit),
    )

    assert model.cell.unit_system is unit_system


@pytest.mark.parametrize("scale", [1.0e-10, 1.0, 1.0e10])
def test__init__rejects_a_nonsquare_cell_independently_of_scale(
    scale: float,
) -> None:
    cell = _cell(
        basis=scale * np.diag([1.0, 1.2]),
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


def test__init__rejects_a_non_bloch_cell() -> None:
    with pytest.raises(TypeError, match="cell must be BlochCell2D"):
        CheckerboardPotentialBloch2D(
            cell=object(),  # type: ignore[arg-type]
            amplitude=ScalarQuantity(1.0, Unitless()),
        )


def test__init__rejects_a_non_scalar_amplitude() -> None:
    with pytest.raises(TypeError, match="amplitude must be ScalarQuantity"):
        CheckerboardPotentialBloch2D(
            cell=_cell(basis=np.eye(2), unit_system=UnitSystem.NONDIMENSIONAL),
            amplitude=1.0,  # type: ignore[arg-type]
        )


def test__init__rejects_an_amplitude_outside_the_selected_system() -> None:
    with pytest.raises(
        ValueError,
        match="amplitude must use the selected energy unit",
    ):
        CheckerboardPotentialBloch2D(
            cell=_cell(basis=np.eye(2), unit_system=UnitSystem.METAL),
            amplitude=ScalarQuantity(1.0, Unitless()),
        )
