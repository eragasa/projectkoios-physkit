import numpy as np
import pytest
from grid_1d import Grid1D, GridType1D


def test_basic_properties():
    a = 0.0
    b = 1.0
    N = 10
    L = b - a
    g = Grid1D(a=a, b=b, N=N, grid_type=GridType1D.LEFT_CLOSED)
    assert g.a == a
    assert g.b == b
    assert g.L == L
    assert g.N == N


def test_closed_grid():
    a = 0
    b = 1
    N = 5
    grid_type = GridType1D.CLOSED

    g = Grid1D(a=a, b=b, N=N, grid_type=grid_type)

    assert isinstance(g.X, np.ndarray)
    assert g.X.shape == (N,)

    assert np.isclose(g.X[0], a)
    assert np.isclose(g.X[-1], b)

    dX = np.diff(g.X)
    assert np.all(dX > 0.0)
    assert np.allclose(dX, g.dx)

    assert np.isclose(g.dx, (b - a) / (N - 1))


def test_closed_grid_requires_N_ge_2():
    with pytest.raises(ValueError):
        Grid1D(a=0.0, b=1.0, N=1.0, grid_type=GridType1D.CLOSED)


def test_left_closed_grid():
    a, b, N = -0.5, 0.5, 10
    g = Grid1D(a, b, N, GridType1D.LEFT_CLOSED)

    assert np.isclose(g.X[0], a)
    assert g.X[-1] < b
    assert np.isclose(g.dx, (b - a) / N)


def test_right_closed_grid():
    a, b, N = -0.5, 0.5, 10
    g = Grid1D(a, b, N, GridType1D.RIGHT_CLOSED)

    assert g.X[0] > a
    assert np.isclose(g.X[-1], b)
    assert np.isclose(g.dx, (b - a) / N)


def test_open_grid():
    a, b, N = -0.5, 0.5, 10
    g = Grid1D(a, b, N, GridType1D.OPEN)

    assert g.X[0] > a
    assert g.X[-1] < b
    assert np.isclose(g.dx, (b - a) / (N + 1))


def test_interior_grid():
    a, b, N = -0.5, 0.5, 10
    g = Grid1D(a, b, N, GridType1D.INTERIOR)

    assert g.X[0] > a
    assert g.X[-1] < b
    assert np.isclose(g.dx, (b - a) / (N + 1))


def test_midpoint_grid():
    a, b, N = -0.5, 0.5, 10
    g = Grid1D(a, b, N, GridType1D.MIDPOINT)

    assert np.isclose(g.X[0], a + 0.5 * g.dx)
    assert np.isclose(g.X[-1], b - 0.5 * g.dx)
    assert np.isclose(g.dx, (b - a) / N)
