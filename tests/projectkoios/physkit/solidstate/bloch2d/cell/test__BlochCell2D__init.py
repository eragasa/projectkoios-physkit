"""Bloch direct-cell construction requirements."""

import numpy as np
import pytest

from projectkoios.physkit.solidstate.bloch2d.cell import BlochCell2D
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
    UnitSystem,
)


def test__init__accepts_an_explicitly_unitless_cell() -> None:
    cell = BlochCell2D(
        direct_basis=MatrixQuantity(
            magnitude=np.array([[1.0, 0.25], [0.0, 0.75]]),
            unit=Unitless(),
        ),
        unit_system=UnitSystem.NONDIMENSIONAL,
    )

    assert isinstance(cell.area.unit, Unitless)
    assert isinstance(cell.reciprocal_basis.unit, Unitless)
    assert isinstance(cell.metric_tensor.unit, Unitless)


def test__init__rejects_a_nonquantity_direct_basis() -> None:
    with pytest.raises(TypeError, match="direct_basis must be MatrixQuantity"):
        BlochCell2D(
            direct_basis=np.eye(2),  # type: ignore[arg-type]
            unit_system=UnitSystem.METAL,
        )


def test__init__rejects_a_non_unit_system() -> None:
    with pytest.raises(TypeError, match="unit_system must be UnitSystem"):
        BlochCell2D(
            direct_basis=MatrixQuantity(
                magnitude=np.eye(2),
                unit=UnitSystem.METAL.length_unit,
            ),
            unit_system="metal",  # type: ignore[arg-type]
        )


def test__init__rejects_the_si_unit_system() -> None:
    with pytest.raises(
        ValueError,
        match="unit_system must be METAL or NONDIMENSIONAL",
    ):
        BlochCell2D(
            direct_basis=MatrixQuantity(
                magnitude=np.eye(2),
                unit=UnitSystem.SI.length_unit,
            ),
            unit_system=UnitSystem.SI,
        )


def test__init__rejects_a_basis_unit_outside_the_selected_system() -> None:
    with pytest.raises(
        ValueError,
        match="direct_basis must use the selected unit-system length unit",
    ):
        BlochCell2D(
            direct_basis=MatrixQuantity(
                magnitude=np.eye(2),
                unit=PhysicalUnit("nanometer"),
            ),
            unit_system=UnitSystem.METAL,
        )


def test__init__rejects_a_basis_with_the_wrong_shape() -> None:
    with pytest.raises(ValueError, match=r"direct_basis must have shape \(2, 2\)"):
        BlochCell2D(
            direct_basis=MatrixQuantity(
                magnitude=np.ones((2, 3)),
                unit=UnitSystem.METAL.length_unit,
            ),
            unit_system=UnitSystem.METAL,
        )


@pytest.mark.parametrize("scale", [1.0e-10, 1.0, 1.0e10])
def test__init__rejects_numerical_rank_deficiency_independently_of_scale(
    scale: float,
) -> None:
    nearly_parallel_basis = scale * np.array(
        [
            [1.0, 1.0],
            [0.0, 0.5e-8],
        ]
    )

    with pytest.raises(
        ValueError,
        match="direct_basis columns must be numerically linearly independent",
    ):
        BlochCell2D(
            direct_basis=MatrixQuantity(
                magnitude=nearly_parallel_basis,
                unit=UnitSystem.METAL.length_unit,
            ),
            unit_system=UnitSystem.METAL,
        )
