"""Verification of ``DirichletBoundaryCondition`` construction."""

import pytest

from projectkoios.physkit.core.boundaries import (
    BoundaryCondition,
    BoundaryConditionType,
)
from projectkoios.physkit.discretization.boundary_conditions import (
    DirichletBoundaryCondition,
)
from projectkoios.physkit.units import ScalarQuantity, Unitless


def test_uses_nominal_boundary_ownership() -> None:
    boundary = DirichletBoundaryCondition(ScalarQuantity(0.0, Unitless()))

    assert isinstance(boundary, BoundaryCondition)
    assert boundary.type is BoundaryConditionType.DIRICHLET
    assert boundary.condition_kind == "dirichlet"


def test_requires_scalar_quantity() -> None:
    with pytest.raises(TypeError):
        DirichletBoundaryCondition(0.0)  # type: ignore[arg-type]
