"""Verification of ``DirichletInterval`` construction."""

import pytest

from projectkoios.physkit.discretization.boundary_conditions import (
    DirichletBoundaryCondition,
)
from projectkoios.physkit.discretization.cartesian_grids import UniformCartesianGrid1D
from projectkoios.physkit.discretization.intervals import DirichletInterval
from projectkoios.physkit.units import ScalarQuantity, Unitless


def grid(spacing: float = 0.5) -> UniformCartesianGrid1D:
    return UniformCartesianGrid1D(
        ScalarQuantity(-1.0, Unitless()),
        ScalarQuantity(1.0, Unitless()),
        ScalarQuantity(spacing, Unitless()),
    )


def test_retains_exact_nominal_representations() -> None:
    represented_grid = grid()
    boundary = DirichletBoundaryCondition(ScalarQuantity(0.0, Unitless()))
    interval = DirichletInterval(represented_grid, boundary)

    assert interval.grid is represented_grid
    assert interval.boundary_condition is boundary
    assert interval.interior_points == 3


def test_rejects_wrong_types_and_empty_interior() -> None:
    boundary = DirichletBoundaryCondition(ScalarQuantity(0.0, Unitless()))
    with pytest.raises(TypeError):
        DirichletInterval(object(), boundary)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        DirichletInterval(grid(), object())  # type: ignore[arg-type]
    empty = UniformCartesianGrid1D(
        ScalarQuantity(0.0, Unitless()),
        ScalarQuantity(1.0, Unitless()),
        ScalarQuantity(1.0, Unitless()),
    )
    with pytest.raises(ValueError, match="interior point"):
        DirichletInterval(empty, boundary)
