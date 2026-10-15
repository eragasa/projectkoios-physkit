"""Software verification for the accepted geometry-only uniform grid."""

from numbers import Integral, Real

import numpy as np
import pytest

from projectkoios.physkit.discretization import UniformGrid1D
from projectkoios.physkit.discretization.grid_1d import (
    UniformGrid1D as DefiningUniformGrid1D,
)


def test_public_import_is_the_defining_class() -> None:
    assert UniformGrid1D is DefiningUniformGrid1D


@pytest.mark.parametrize(
    "a,b,num_points",
    [
        (0, 1, 3),
        (np.float32(-2.5), np.float64(1.5), np.int64(9)),
    ],
)
def test_valid_scalar_inputs_are_canonicalized_and_geometry_is_correct(
    a: Real, b: Real, num_points: Integral
) -> None:
    grid = UniformGrid1D(a, b, num_points)

    assert type(grid.a) is float
    assert type(grid.b) is float
    assert type(grid.num_points) is int
    assert grid.length == float(b) - float(a)
    assert grid.spacing == (float(b) - float(a)) / (int(num_points) - 1)
    np.testing.assert_array_equal(
        grid.coordinates,
        np.linspace(float(a), float(b), int(num_points), dtype=np.float64),
    )


@pytest.mark.parametrize(
    "name", ["boundary_values", "active_indices", "restrict", "embed"]
)
def test_grid_exposes_geometry_only(name: str) -> None:
    assert not hasattr(UniformGrid1D(0.0, 1.0, 5), name)


def test_coordinates_are_owned_c_contiguous_float64_copies() -> None:
    grid = UniformGrid1D(-1.0, 2.0, 7)
    first = grid.coordinates
    second = grid.coordinates

    assert first.dtype == np.dtype(np.float64)
    assert first.flags.c_contiguous
    assert first.flags.owndata
    assert not np.shares_memory(first, second)
    first[:] = 99.0
    np.testing.assert_array_equal(
        grid.coordinates, np.linspace(-1.0, 2.0, 7, dtype=np.float64)
    )


@pytest.mark.parametrize(
    "attribute,value", [("a", 2.0), ("b", 3.0), ("num_points", 8), ("extra", 1)]
)
def test_grid_is_observably_immutable(attribute: str, value: object) -> None:
    grid = UniformGrid1D(0.0, 1.0, 5)

    with pytest.raises(AttributeError):
        setattr(grid, attribute, value)


@pytest.mark.parametrize(
    "a,b,num_points",
    [
        (True, 1.0, 3),
        (np.bool_(False), 1.0, 3),
        ("0", 1.0, 3),
        (0.0, False, 3),
        (0.0, 1.0 + 0.0j, 3),
        (0.0, 1.0, True),
        (0.0, 1.0, np.float64(3.0)),
        (0.0, 1.0, "3"),
    ],
)
def test_wrong_scalar_semantic_types_raise_type_error(
    a: object, b: object, num_points: object
) -> None:
    with pytest.raises(TypeError):
        UniformGrid1D(a, b, num_points)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "a,b,num_points",
    [
        (np.nan, 1.0, 3),
        (0.0, np.inf, 3),
        (1.0, 1.0, 3),
        (2.0, 1.0, 3),
        (0.0, 1.0, 2),
        (0.0, 1.0, -1),
        (10**400, 1.0, 3),
        (-(10**400), 1.0, 3),
        (0, 10**400, 3),
        (0, -(10**400), 3),
    ],
)
def test_correct_scalar_types_violating_invariants_raise_value_error(
    a: Real, b: Real, num_points: Integral
) -> None:
    with pytest.raises(ValueError):
        UniformGrid1D(a, b, num_points)
