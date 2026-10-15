"""Software and numerical verification for the accepted 1D Laplacian."""

from typing import get_type_hints

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from projectkoios.physkit.discretization import (
    HomogeneousDirichletStateSpace1D,
    UniformGrid1D,
)
from projectkoios.physkit.operators import (
    DiscreteLinearOperator1D,
    FiniteDifferenceLaplacian1D,
    LinearOperator,
)
from projectkoios.physkit.operators.finite_difference_1d import (
    FiniteDifferenceLaplacian1D as DefiningFiniteDifferenceLaplacian1D,
)

RTOL = 5e-14


def make_operator(
    num_points: int, a: float = 0.0, b: float = 1.0
) -> FiniteDifferenceLaplacian1D:
    grid = UniformGrid1D(a, b, num_points)
    return FiniteDifferenceLaplacian1D(HomogeneousDirichletStateSpace1D(grid))


def direct_tridiagonal_oracle(num_points: int, a: float, b: float) -> np.ndarray:
    """Construct the accepted matrix entry-by-entry without sparse helpers."""
    dimension = num_points - 2
    spacing = (b - a) / (num_points - 1)
    expected = np.zeros((dimension, dimension), dtype=np.float64)
    for row in range(dimension):
        expected[row, row] = -2.0 / spacing**2
        if row > 0:
            expected[row, row - 1] = 1.0 / spacing**2
        if row + 1 < dimension:
            expected[row, row + 1] = 1.0 / spacing**2
    return expected


def componentwise_stencil_oracle(state: np.ndarray, spacing: float) -> np.ndarray:
    """Apply the centered stencil to explicitly zero-padded endpoint data."""
    full = np.zeros(state.size + 2, dtype=state.dtype)
    full[1:-1] = state
    expected = np.empty(state.size, dtype=state.dtype)
    for active_index in range(1, full.size - 1):
        expected[active_index - 1] = (
            full[active_index - 1] - 2.0 * full[active_index] + full[active_index + 1]
        ) / spacing**2
    return expected


def accepted_allclose(actual: np.ndarray, reference: np.ndarray) -> None:
    atol = RTOL * max(1.0, float(np.linalg.norm(reference, ord=np.inf)))
    np.testing.assert_allclose(actual, reference, rtol=RTOL, atol=atol)


def test_public_import_is_defining_class_and_hierarchy_is_exact() -> None:
    assert FiniteDifferenceLaplacian1D is DefiningFiniteDifferenceLaplacian1D
    assert issubclass(FiniteDifferenceLaplacian1D, DiscreteLinearOperator1D)
    assert issubclass(FiniteDifferenceLaplacian1D, LinearOperator)


def test_constructor_requires_homogeneous_dirichlet_state_space() -> None:
    for value in (None, UniformGrid1D(0.0, 1.0, 5), object()):
        with pytest.raises(TypeError):
            FiniteDifferenceLaplacian1D(value)  # type: ignore[arg-type]


def test_domain_codomain_shape_dtype_and_narrow_public_annotations() -> None:
    operator = make_operator(8, -2.0, 3.0)
    domain_hints = get_type_hints(FiniteDifferenceLaplacian1D.domain.fget)  # type: ignore[arg-type]
    codomain_hints = get_type_hints(FiniteDifferenceLaplacian1D.codomain.fget)  # type: ignore[arg-type]

    assert operator.domain is operator.codomain
    assert operator.domain.semantic_identity == operator.codomain.semantic_identity
    assert operator.shape == (6, 6)
    assert operator.dtype == np.dtype(np.float64)
    assert domain_hints["return"] is HomogeneousDirichletStateSpace1D
    assert codomain_hints["return"] is HomogeneousDirichletStateSpace1D


@pytest.mark.parametrize("num_points", [3, 4, 8])
def test_csr_and_dense_entries_equal_independent_direct_oracle(num_points: int) -> None:
    a, b = -0.75, 1.25
    operator = make_operator(num_points, a, b)
    sparse = operator.matrix
    expected = direct_tridiagonal_oracle(num_points, a, b)

    assert isinstance(sparse, csr_matrix)
    assert sparse.format == "csr"
    assert sparse.shape == (num_points - 2, num_points - 2)
    assert sparse.dtype == np.dtype(np.float64)
    assert sparse.has_canonical_format
    assert sparse.has_sorted_indices
    assert sparse.nnz == np.count_nonzero(expected)
    assert np.all(sparse.data != 0.0)
    np.testing.assert_array_equal(sparse.toarray(), expected)
    np.testing.assert_array_equal(operator.to_dense(), expected)


def test_n_equals_three_is_exact_one_by_one_positive_second_derivative_matrix() -> None:
    operator = make_operator(3, 0.0, 1.0)

    np.testing.assert_array_equal(operator.to_dense(), np.array([[-8.0]]))


def test_matrix_is_symmetric_negative_definite_with_positive_off_diagonals() -> None:
    dense = make_operator(8).to_dense()

    np.testing.assert_array_equal(dense, dense.T)
    assert np.all(np.diag(dense) < 0.0)
    assert np.all(np.diag(dense, 1) > 0.0)
    assert np.all(np.linalg.eigvalsh(dense) < 0.0)


def test_sparse_and_dense_inspection_are_defensive_owned_copies() -> None:
    operator = make_operator(8)
    first_sparse = operator.matrix
    second_sparse = operator.matrix
    dense = operator.to_dense()
    expected = direct_tridiagonal_oracle(8, 0.0, 1.0)

    assert not np.shares_memory(first_sparse.data, second_sparse.data)
    first_sparse.data[:] = 0.0
    dense[:] = 0.0
    np.testing.assert_array_equal(operator.matrix.toarray(), expected)
    np.testing.assert_array_equal(operator.to_dense(), expected)
    assert operator.to_dense().flags.c_contiguous
    assert operator.to_dense().flags.owndata


@pytest.mark.parametrize(
    "state",
    [
        np.array([1.0, -2.0, 0.5, 4.0, -3.0, 2.0]),
        np.array([1.0 + 2.0j, -2.0j, 0.5 - 0.25j, 4.0j, -3.0, 2.0 + 7.0j]),
    ],
)
def test_real_and_complex_application_equal_independent_stencil(
    state: np.ndarray,
) -> None:
    operator = make_operator(8, -1.0, 2.0)
    before = state.copy()
    expected = componentwise_stencil_oracle(state, operator.domain.grid.spacing)
    applied = operator.apply(state)
    multiplied = operator @ state

    accepted_allclose(applied, expected)
    accepted_allclose(multiplied, expected)
    np.testing.assert_array_equal(state, before)
    expected_dtype = np.dtype(np.complex128 if np.iscomplexobj(state) else np.float64)
    assert applied.dtype == expected_dtype
    assert applied.flags.c_contiguous and applied.flags.owndata
    assert not np.shares_memory(applied, state)


def test_complex_application_equals_separate_real_and_imaginary_actions() -> None:
    operator = make_operator(8)
    state = np.array([1 + 2j, 2 - 3j, -4 + 0.5j, 5j, 6 - 2j, -1j])

    combined = operator @ state
    separate = (operator @ state.real) + 1j * (operator @ state.imag)
    accepted_allclose(combined, separate)


@pytest.mark.parametrize(
    "invalid",
    [
        [True] * 6,
        [1.0, 2.0, 3.0, False, 5.0, 6.0],
        np.array([1, 2, 3, 4, 5, 6], dtype=object),
        ["1", "2", "3", "4", "5", "6"],
        [[1.0], [2.0, 3.0]],
    ],
)
def test_application_rejects_wrong_semantic_vector_types(invalid: object) -> None:
    operator = make_operator(8)
    with pytest.raises(TypeError):
        operator.apply(invalid)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        operator @ invalid  # type: ignore[operator]


@pytest.mark.parametrize(
    "invalid",
    [np.zeros((6, 1)), np.zeros(5), [1.0, 2.0, 3.0, np.nan, 5.0, 6.0]],
)
def test_application_rank_shape_and_finiteness_violations_raise_value_error(
    invalid: object,
) -> None:
    operator = make_operator(8)
    with pytest.raises(ValueError):
        operator.apply(invalid)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        operator @ invalid  # type: ignore[operator]


@pytest.mark.parametrize("factor", [-0.5, 1.0 + 2.0j])
@pytest.mark.parametrize(
    "state",
    [
        np.array([1.0, -2.0, 0.5, 4.0, -3.0, 2.0]),
        np.array([1.0 + 2.0j, -2.0j, 0.5 - 0.25j, 4.0j, -3.0, 2.0 + 7.0j]),
    ],
)
def test_real_and_complex_scaling_matches_explicit_multiplication(
    factor: complex, state: np.ndarray
) -> None:
    operator = make_operator(8)
    scaled = operator.scaled(factor)
    reference = factor * componentwise_stencil_oracle(
        state, operator.domain.grid.spacing
    )
    result = scaled @ state

    assert scaled.domain is operator.domain
    assert scaled.codomain is operator.codomain
    assert scaled.shape == operator.shape
    assert scaled.dtype == np.dtype(np.result_type(operator.dtype, factor))
    accepted_allclose(result, reference)
    assert result.flags.c_contiguous and result.flags.owndata
    assert not np.shares_memory(result, state)


@pytest.mark.parametrize("mode", [1, 2])
def test_sine_modes_have_monotone_error_and_accepted_final_orders(mode: int) -> None:
    point_counts = (17, 33, 65, 129)
    errors: list[float] = []
    spacings: list[float] = []

    for num_points in point_counts:
        operator = make_operator(num_points)
        x_active = operator.domain.active_coordinates
        spacing = operator.domain.grid.spacing
        sampled = np.sin(mode * np.pi * x_active)
        exact_second_derivative = -((mode * np.pi) ** 2) * sampled
        residual = operator @ sampled - exact_second_derivative
        weighted_residual_norm = np.sqrt(spacing * np.sum(np.abs(residual) ** 2))
        weighted_reference_norm = np.sqrt(
            spacing * np.sum(np.abs(exact_second_derivative) ** 2)
        )
        errors.append(float(weighted_residual_norm / weighted_reference_norm))
        spacings.append(spacing)

    assert all(fine < coarse for coarse, fine in zip(errors, errors[1:], strict=False))
    for coarse_index, fine_index in ((1, 2), (2, 3)):
        order = np.log(errors[coarse_index] / errors[fine_index]) / np.log(
            spacings[coarse_index] / spacings[fine_index]
        )
        assert 1.90 <= order <= 2.10


def test_operator_is_observably_immutable_and_exposes_no_composition() -> None:
    operator = make_operator(8)

    for attribute, value in (("domain", None), ("shape", (1, 1)), ("extra", 1)):
        with pytest.raises(AttributeError):
            setattr(operator, attribute, value)
    assert not hasattr(operator, "compose")
