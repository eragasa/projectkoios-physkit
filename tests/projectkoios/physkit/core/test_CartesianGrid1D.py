import numpy as np
import pytest

from projectkoios.physkit.core.grids import CartesianAxis, CartesianGrid1D


def test__init__stores_axis_object():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
        endpoint=True,
    )

    assert isinstance(grid.x_axis, CartesianAxis)


def test__init__stores_constructor_arguments():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
        endpoint=True,
    )

    assert grid.x_min == 0.0
    assert grid.x_max == 1.0
    assert grid.Nx == 5
    assert grid.endpoint is True


def test__init__builds_x_axis_values_with_endpoint():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
        endpoint=True,
    )

    expected = np.array([0.0, 0.25, 0.50, 0.75, 1.0])

    np.testing.assert_allclose(grid.x, expected)


def test__init__builds_x_axis_values_without_endpoint():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
        endpoint=False,
    )

    expected = np.array([0.0, 0.2, 0.4, 0.6, 0.8])

    np.testing.assert_allclose(grid.x, expected)


def test__init__computes_domain_length():
    grid = CartesianGrid1D(
        x_min=-2.0,
        x_max=3.0,
        Nx=6,
        endpoint=True,
    )

    assert grid.Lx == pytest.approx(5.0)


def test__init__computes_spacing_with_endpoint():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
        endpoint=True,
    )

    assert grid.dx == pytest.approx(0.25)


def test__init__computes_spacing_without_endpoint():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
        endpoint=False,
    )

    assert grid.dx == pytest.approx(0.2)


def test__init__sets_shape():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
    )

    assert grid.shape == (5,)


def test__init__sets_size():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
    )

    assert grid.size == 5


def test__init__x_is_float64():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
    )

    assert grid.x.dtype == np.float64


def test__init__x_is_read_only():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
    )

    assert grid.x.flags.writeable is False

    with pytest.raises(ValueError):
        grid.x[0] = 10.0


def test__init__x_references_axis_values():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
    )

    assert grid.x is grid.x_axis.values


@pytest.mark.parametrize("bad_x_min", [np.nan, np.inf, -np.inf])
def test__init__rejects_nonfinite_x_min(bad_x_min):
    with pytest.raises(ValueError, match="x_min must be finite"):
        CartesianGrid1D(
            x_min=bad_x_min,
            x_max=1.0,
            Nx=5,
        )


@pytest.mark.parametrize("bad_x_max", [np.nan, np.inf, -np.inf])
def test__init__rejects_nonfinite_x_max(bad_x_max):
    with pytest.raises(ValueError, match="x_max must be finite"):
        CartesianGrid1D(
            x_min=0.0,
            x_max=bad_x_max,
            Nx=5,
        )


@pytest.mark.parametrize(
    "x_min, x_max",
    [
        (1.0, 1.0),
        (2.0, 1.0),
    ],
)
def test__init__rejects_invalid_bounds(x_min, x_max):
    with pytest.raises(ValueError, match="x_max must be greater than x_min"):
        CartesianGrid1D(
            x_min=x_min,
            x_max=x_max,
            Nx=5,
        )


@pytest.mark.parametrize("bad_Nx", [5.0, "5", None])
def test__init__rejects_noninteger_Nx(bad_Nx):
    with pytest.raises(TypeError, match="Nx must be an integer"):
        CartesianGrid1D(
            x_min=0.0,
            x_max=1.0,
            Nx=bad_Nx,
        )


@pytest.mark.parametrize("bad_Nx", [-1, 0, 1])
def test__init__rejects_too_few_points(bad_Nx):
    with pytest.raises(ValueError, match="Nx must be at least 2"):
        CartesianGrid1D(
            x_min=0.0,
            x_max=1.0,
            Nx=bad_Nx,
        )


@pytest.mark.parametrize("bad_endpoint", [0, 1, "True", None])
def test__init__rejects_nonbool_endpoint(bad_endpoint):
    with pytest.raises(TypeError, match="endpoint_x must be a bool"):
        CartesianGrid1D(
            x_min=0.0,
            x_max=1.0,
            Nx=5,
            endpoint=bad_endpoint,
        )


def test__repr__returns_expected_string():
    grid = CartesianGrid1D(
        x_min=0.0,
        x_max=1.0,
        Nx=5,
        endpoint=True,
    )

    expected = "CartesianGrid1D(x_min=0.0, x_max=1.0, Nx=5, endpoint=True)"

    assert repr(grid) == expected
