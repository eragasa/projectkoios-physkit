"""Construction tests for the bounded oblique PIAB2D model."""

import numpy as np
import pytest

from projectkoios.physkit.qm.piab2d.oblique import ObliquePiab2D
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    UnitSystem,
)


def test_accepts_unitless_linearly_independent_direct_basis() -> None:
    """A finite skew basis defines a positive-area bounded model."""
    model = ObliquePiab2D(
        direct_basis=MatrixQuantity(
            np.array([[1.0, 0.4], [0.0, 0.9]]),
            Unitless(),
        ),
        mass=ScalarQuantity(1.0, Unitless()),
        unit_system=UnitSystem.NONDIMENSIONAL,
    )

    assert model.dimension == 2
    assert model.area.magnitude == pytest.approx(0.9)
    assert isinstance(model.area.unit, Unitless)
    assert model.metric.inverse_metric_tensor[0, 1] != 0.0


@pytest.mark.parametrize(
    ("basis", "match"),
    [
        (np.ones((2, 3)), "shape"),
        (np.array([[1.0, 2.0], [0.0, 0.0]]), "linearly independent"),
    ],
)
def test_rejects_invalid_direct_geometry(
    basis: np.ndarray,
    match: str,
) -> None:
    """Wrong shape and collinear edge vectors cannot define a cell."""
    with pytest.raises(ValueError, match=match):
        ObliquePiab2D(
            direct_basis=MatrixQuantity(basis, Unitless()),
            mass=ScalarQuantity(1.0, Unitless()),
            unit_system=UnitSystem.NONDIMENSIONAL,
        )


def test_rejects_unit_system_mismatches() -> None:
    """Length and mass declarations must use one coherent selected system."""
    with pytest.raises(ValueError, match="selected length unit"):
        ObliquePiab2D(
            direct_basis=MatrixQuantity(np.eye(2), PhysicalUnit("meter")),
            mass=ScalarQuantity(1.0, Unitless()),
            unit_system=UnitSystem.NONDIMENSIONAL,
        )
    with pytest.raises(ValueError, match="selected mass unit"):
        ObliquePiab2D(
            direct_basis=MatrixQuantity(np.eye(2), Unitless()),
            mass=ScalarQuantity(1.0, PhysicalUnit("kilogram")),
            unit_system=UnitSystem.NONDIMENSIONAL,
        )


def test_rejects_nonpositive_mass() -> None:
    """A free-particle kinetic model requires positive mass."""
    with pytest.raises(ValueError, match="positive"):
        ObliquePiab2D(
            direct_basis=MatrixQuantity(np.eye(2), Unitless()),
            mass=ScalarQuantity(0.0, Unitless()),
            unit_system=UnitSystem.NONDIMENSIONAL,
        )
