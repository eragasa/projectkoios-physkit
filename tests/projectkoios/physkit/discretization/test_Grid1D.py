import numpy as np
import pytest

from projectkoios.physkit.discretization.grid_1d import ActiveSetType1D, Grid1D


def test__Grid1D__x_returns_closed_uniform_grid():
    grid = Grid1D(a=0.0, b=1.0, N=5)

    expected = np.linspace(0.0, 1.0, 5)

    assert np.allclose(grid.x, expected)


def test__Grid1D__L_returns_domain_length():
    grid = Grid1D(a=-1.0, b=2.0, N=5)

    assert np.isclose(grid.L, 3.0)


def test__Grid1D__dx_returns_closed_grid_spacing():
    grid = Grid1D(a=0.0, b=1.0, N=5)

    assert np.isclose(grid.dx, 0.25)


@pytest.mark.parametrize(
    "active_type, expected",
    [
        (ActiveSetType1D.ALL, np.arange(0, 10)),
        (ActiveSetType1D.INTERIOR, np.arange(1, 9)),
        (ActiveSetType1D.LEFT_CLOSED, np.arange(0, 9)),
        (ActiveSetType1D.RIGHT_CLOSED, np.arange(1, 10)),
        (ActiveSetType1D.LEFT_BOUNDARY, np.array([0])),
        (ActiveSetType1D.RIGHT_BOUNDARY, np.array([9])),
        (ActiveSetType1D.BOUNDARY, np.array([0, 9])),
    ],
)
def test__Grid1D__active_indices_match_active_set_type(active_type, expected):
    grid = Grid1D(a=0.0, b=1.0, N=10, active_type=active_type)

    assert np.array_equal(grid.active_indices, expected)


def test__Grid1D__x_active_returns_active_coordinates():
    grid = Grid1D(
        a=0.0,
        b=1.0,
        N=10,
        active_type=ActiveSetType1D.INTERIOR,
    )

    expected = grid.x[1:9]

    assert np.allclose(grid.x_active, expected)


def test__Grid1D__constructor_raises_for_invalid_domain():
    with pytest.raises(ValueError, match="Require b > a"):
        Grid1D(a=1.0, b=0.0, N=10)


def test__Grid1D__constructor_raises_for_too_few_points():
    with pytest.raises(ValueError, match="Require N >= 2"):
        Grid1D(a=0.0, b=1.0, N=1)


def test__Grid1D__interior_active_indices_empty_for_two_point_grid():
    grid = Grid1D(
        a=0.0,
        b=1.0,
        N=2,
        active_type=ActiveSetType1D.INTERIOR,
    )

    assert np.array_equal(grid.active_indices, np.array([]))


def test__Grid1D__boundary_active_indices_for_two_point_grid():
    grid = Grid1D(
        a=0.0,
        b=1.0,
        N=2,
        active_type=ActiveSetType1D.BOUNDARY,
    )

    assert np.array_equal(grid.active_indices, np.array([0, 1]))
