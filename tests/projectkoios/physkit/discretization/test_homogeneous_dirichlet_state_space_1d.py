"""Software verification for homogeneous-Dirichlet state-space semantics."""

import numpy as np
import pytest

from projectkoios.physkit.discretization import (
    HomogeneousDirichletStateSpace1D,
    UniformGrid1D,
)
from projectkoios.physkit.discretization.state_space_1d import (
    HomogeneousDirichletStateSpace1D as DefiningStateSpace,
)


def make_space(num_points: int = 5) -> HomogeneousDirichletStateSpace1D:
    return HomogeneousDirichletStateSpace1D(UniformGrid1D(-1.0, 2.0, num_points))


def test_public_import_is_the_defining_class() -> None:
    assert HomogeneousDirichletStateSpace1D is DefiningStateSpace


def test_constructor_requires_uniform_grid() -> None:
    for value in (None, object(), [0.0, 0.5, 1.0]):
        with pytest.raises(TypeError):
            HomogeneousDirichletStateSpace1D(value)  # type: ignore[arg-type]


def test_state_space_exposes_exact_geometry_boundary_and_active_data() -> None:
    grid = UniformGrid1D(-1.0, 2.0, 5)
    space = HomogeneousDirichletStateSpace1D(grid)

    assert space.grid is grid
    assert space.boundary_values == (0.0, 0.0)
    assert space.dimension == 3
    np.testing.assert_array_equal(
        space.active_indices, np.array([1, 2, 3], dtype=np.intp)
    )
    np.testing.assert_array_equal(space.active_coordinates, grid.coordinates[1:-1])
    assert space.active_indices.dtype == np.dtype(np.intp)
    assert space.active_coordinates.dtype == np.dtype(np.float64)


def test_semantic_identity_is_exact_and_value_based() -> None:
    first = make_space()
    equivalent = make_space()
    different = HomogeneousDirichletStateSpace1D(UniformGrid1D(-1.0, 2.0, 6))
    expected = (
        "HomogeneousDirichletStateSpace1D",
        ("UniformGrid1D", -1.0, 2.0, 5, "closed-endpoint-inclusive"),
        "active-interior-indices-1-through-N-minus-2",
        ("homogeneous-dirichlet", 0.0, 0.0),
        "real-or-complex-active-vectors",
    )

    assert first.semantic_identity == expected
    assert first.semantic_identity == equivalent.semantic_identity
    assert first.semantic_identity != different.semantic_identity
    assert first is not equivalent


@pytest.mark.parametrize("property_name", ["active_indices", "active_coordinates"])
def test_array_properties_return_owned_defensive_copies(property_name: str) -> None:
    space = make_space()
    first = getattr(space, property_name)
    second = getattr(space, property_name)

    assert first.flags.c_contiguous
    assert first.flags.owndata
    assert not np.shares_memory(first, second)
    first[:] = 0
    assert not np.array_equal(first, getattr(space, property_name))


@pytest.mark.parametrize(
    "attribute,value", [("grid", None), ("dimension", 99), ("extra", 1)]
)
def test_state_space_is_observably_immutable(attribute: str, value: object) -> None:
    space = make_space()

    with pytest.raises(AttributeError):
        setattr(space, attribute, value)


def test_restrict_discards_unchecked_endpoints_and_owns_real_result() -> None:
    space = make_space()
    full = np.array([123.0, 1.5, -2.0, 4.25, -987.0], dtype=np.float32)
    before = full.copy()
    restricted = space.restrict(full)

    np.testing.assert_array_equal(restricted, np.array([1.5, -2.0, 4.25]))
    np.testing.assert_array_equal(full, before)
    assert restricted.dtype == np.dtype(np.float64)
    assert restricted.flags.c_contiguous
    assert restricted.flags.owndata
    assert not np.shares_memory(restricted, full)


def test_restrict_and_embed_preserve_complex_values_in_owned_outputs() -> None:
    space = make_space()
    full = np.array([9 + 2j, 1 + 3j, -2j, 4 - 5j, 7j], dtype=np.complex64)
    restricted = space.restrict(full)
    embedded = space.embed(restricted)

    np.testing.assert_array_equal(restricted, full[1:-1].astype(np.complex128))
    np.testing.assert_array_equal(
        embedded, np.array([0, 1 + 3j, -2j, 4 - 5j, 0], dtype=np.complex128)
    )
    assert restricted.dtype == np.dtype(np.complex128)
    assert embedded.dtype == np.dtype(np.complex128)
    assert restricted.flags.c_contiguous and restricted.flags.owndata
    assert embedded.flags.c_contiguous and embedded.flags.owndata
    assert embedded[0] == 0j and embedded[-1] == 0j


def test_embed_owns_real_result_and_does_not_mutate_or_retain_input() -> None:
    space = make_space()
    active = np.array([1, 2, 3], dtype=np.int16)
    embedded = space.embed(active)
    active[:] = 9

    np.testing.assert_array_equal(embedded, np.array([0.0, 1.0, 2.0, 3.0, 0.0]))
    assert embedded.dtype == np.dtype(np.float64)
    assert embedded.flags.c_contiguous and embedded.flags.owndata
    assert not np.shares_memory(embedded, active)


@pytest.mark.parametrize(
    "invalid",
    [
        [True, False, True, False, True],
        [0.0, 1.0, False, 3.0, 4.0],
        np.array([1, 2, 3, 4, 5], dtype=object),
        ["0", "1", "2", "3", "4"],
        [[0.0], [1.0, 2.0]],
        [None, 1.0, 2.0, 3.0, 4.0],
    ],
)
def test_restrict_rejects_wrong_semantic_vector_types(invalid: object) -> None:
    with pytest.raises(TypeError):
        make_space().restrict(invalid)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "invalid",
    [
        [True, False, True],
        [1.0, np.bool_(False), 3.0],
        np.array([1, 2, 3], dtype=object),
        ["1", "2", "3"],
        [[1.0], [2.0, 3.0]],
    ],
)
def test_embed_rejects_wrong_semantic_vector_types(invalid: object) -> None:
    with pytest.raises(TypeError):
        make_space().embed(invalid)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "method,value",
    [
        ("restrict", np.zeros((5, 1))),
        ("restrict", np.zeros(4)),
        ("restrict", [0.0, 1.0, np.nan, 3.0, 4.0]),
        ("restrict", [0.0, 1.0, np.inf, 3.0, 4.0]),
        ("embed", np.zeros((3, 1))),
        ("embed", [[1.0], [2.0], [3.0]]),
        ("embed", np.zeros(2)),
        ("embed", [1.0, np.nan, 3.0]),
        ("embed", [1.0, np.inf, 3.0]),
    ],
)
def test_rank_shape_and_finiteness_violations_raise_value_error(
    method: str, value: object
) -> None:
    with pytest.raises(ValueError):
        getattr(make_space(), method)(value)
