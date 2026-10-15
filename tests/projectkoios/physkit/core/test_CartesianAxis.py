import numpy as np
import pytest

from projectkoios.physkit.core.grids import CartesianAxis


def test__init__stores_constructor_arguments():
    axis = CartesianAxis(
        name="x",
        lower=0.0,
        upper=1.0,
        N=5,
        endpoint=True,
    )

    assert axis.name == "x"
    assert axis.lower == 0.0
    assert axis.upper == 1.0
    assert axis.N == 5
    assert axis.endpoint is True


def test__init__builds_values_with_endpoint():
    axis = CartesianAxis(
        name="x",
        lower=0.0,
        upper=1.0,
        N=5,
        endpoint=True,
    )

    expected = np.array([0.0, 0.25, 0.50, 0.75, 1.0])

    np.testing.assert_allclose(axis.values, expected)


def test__init__builds_values_without_endpoint():
    axis = CartesianAxis(
        name="x",
        lower=0.0,
        upper=1.0,
        N=5,
        endpoint=False,
    )

    expected = np.array([0.0, 0.2, 0.4, 0.6, 0.8])

    np.testing.assert_allclose(axis.values, expected)


def test__init__computes_length():
    axis = CartesianAxis(
        name="x",
        lower=-2.0,
        upper=3.0,
        N=6,
        endpoint=True,
    )

    assert pytest.approx(5.0) == axis.L


def test__init__computes_spacing_with_endpoint():
    axis = CartesianAxis(
        name="x",
        lower=0.0,
        upper=1.0,
        N=5,
        endpoint=True,
    )

    assert axis.d == pytest.approx(0.25)


def test__init__computes_spacing_without_endpoint():
    axis = CartesianAxis(
        name="x",
        lower=0.0,
        upper=1.0,
        N=5,
        endpoint=False,
    )

    assert axis.d == pytest.approx(0.2)


def test__init__makes_values_float64():
    axis = CartesianAxis(
        name="x",
        lower=0.0,
        upper=1.0,
        N=5,
    )

    assert axis.values.dtype == np.float64


def test__init__makes_values_read_only():
    axis = CartesianAxis(
        name="x",
        lower=0.0,
        upper=1.0,
        N=5,
    )

    assert axis.values.flags.writeable is False

    with pytest.raises(ValueError):
        axis.values[0] = 10.0


def test__init__rejects_non_string_name():
    with pytest.raises(TypeError, match="name must be a string"):
        CartesianAxis(
            name=1,
            lower=0.0,
            upper=1.0,
            N=5,
        )


@pytest.mark.parametrize("bad_lower", [np.nan, np.inf, -np.inf])
def test__init__rejects_nonfinite_lower(bad_lower):
    with pytest.raises(ValueError, match="x_min must be finite"):
        CartesianAxis(
            name="x",
            lower=bad_lower,
            upper=1.0,
            N=5,
        )


@pytest.mark.parametrize("bad_upper", [np.nan, np.inf, -np.inf])
def test__init__rejects_nonfinite_upper(bad_upper):
    with pytest.raises(ValueError, match="x_max must be finite"):
        CartesianAxis(
            name="x",
            lower=0.0,
            upper=bad_upper,
            N=5,
        )


@pytest.mark.parametrize(
    "lower, upper",
    [
        (1.0, 1.0),
        (2.0, 1.0),
    ],
)
def test__init__rejects_invalid_bounds(lower, upper):
    with pytest.raises(ValueError, match="x_max must be greater than x_min"):
        CartesianAxis(
            name="x",
            lower=lower,
            upper=upper,
            N=5,
        )


@pytest.mark.parametrize("bad_N", [5.0, "5", None])
def test__init__rejects_noninteger_N(bad_N):
    with pytest.raises(TypeError, match="Nx must be an integer"):
        CartesianAxis(
            name="x",
            lower=0.0,
            upper=1.0,
            N=bad_N,
        )


@pytest.mark.parametrize("bad_N", [-1, 0, 1])
def test__init__rejects_too_few_points(bad_N):
    with pytest.raises(ValueError, match="Nx must be at least 2"):
        CartesianAxis(
            name="x",
            lower=0.0,
            upper=1.0,
            N=bad_N,
        )


@pytest.mark.parametrize("bad_endpoint", [0, 1, "True", None])
def test__init__rejects_nonbool_endpoint(bad_endpoint):
    with pytest.raises(TypeError, match="endpoint_x must be a bool"):
        CartesianAxis(
            name="x",
            lower=0.0,
            upper=1.0,
            N=5,
            endpoint=bad_endpoint,
        )


def test__repr__returns_expected_string():
    axis = CartesianAxis(
        name="x",
        lower=0.0,
        upper=1.0,
        N=5,
        endpoint=True,
    )

    expected = "CartesianAxis(name='x', lower=0.0, upper=1.0, N=5, endpoint=True)"

    assert repr(axis) == expected
