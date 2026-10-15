"""Software verification for the accepted public operator abstractions."""

from typing import get_type_hints

import numpy as np
import pytest
from scipy.sparse import csr_matrix

import projectkoios.physkit.operators as public_operators
from projectkoios.physkit.discretization import (
    HomogeneousDirichletStateSpace1D,
    UniformGrid1D,
)
from projectkoios.physkit.operators import DiscreteLinearOperator1D, LinearOperator
from projectkoios.physkit.operators.base import LinearOperator as DefiningLinearOperator
from projectkoios.physkit.operators.discrete_1d import (
    DiscreteLinearOperator1D as DefiningDiscreteLinearOperator1D,
)


class DiagonalOperator(DiscreteLinearOperator1D):
    """Public-contract-only test double with immutable tuple endpoints."""

    def __init__(self) -> None:
        self._domain = ("domain", 2)
        self._codomain = ("codomain", 2)
        self._matrix = csr_matrix(np.diag([2.0, -3.0]))

    @property
    def domain(self) -> object:
        return self._domain

    @property
    def codomain(self) -> object:
        return self._codomain

    @property
    def shape(self) -> tuple[int, int]:
        return (2, 2)

    @property
    def dtype(self) -> np.dtype:
        return np.dtype(np.float64)

    @property
    def matrix(self) -> csr_matrix:
        return self._matrix.copy()

    def apply(self, state: object) -> np.ndarray:
        vector = np.array(state, dtype=np.float64, copy=True)
        if vector.shape != (2,):
            raise ValueError("state must have shape (2,).")
        return np.array(self._matrix @ vector, order="C", copy=True)


def test_public_imports_are_defining_classes() -> None:
    assert LinearOperator is DefiningLinearOperator
    assert DiscreteLinearOperator1D is DefiningDiscreteLinearOperator1D


def test_public_operator_surface_contains_only_accepted_nonsymbolic_api() -> None:
    expected = {
        "LinearOperator",
        "DiscreteLinearOperator1D",
        "FiniteDifferenceLaplacian1D",
    }
    forbidden = {
        "compose",
        "SymbolicOperator",
        "EquationSpecification",
        "EquationSpec",
        "EquationRegistry",
        "equation_registry",
        "equation_decorator",
        "EquationRenderer",
        "EquationCatalog",
        "QuantumKineticEnergy1D",
        "Hamiltonian1D",
        "Hamiltonian",
    }

    assert set(public_operators.__all__) == expected
    assert forbidden.isdisjoint(public_operators.__all__)
    assert all(not hasattr(public_operators, name) for name in forbidden)


def test_base_classes_are_abstract_and_discrete_inherits_linear() -> None:
    with pytest.raises(TypeError):
        LinearOperator()
    with pytest.raises(TypeError):
        DiscreteLinearOperator1D()
    assert issubclass(DiscreteLinearOperator1D, LinearOperator)


def test_general_domain_and_codomain_annotations_remain_object_only() -> None:
    domain_hints = get_type_hints(LinearOperator.domain.fget)  # type: ignore[arg-type]
    codomain_hints = get_type_hints(LinearOperator.codomain.fget)  # type: ignore[arg-type]

    assert domain_hints["return"] is object
    assert codomain_hints["return"] is object


def test_matmul_delegates_to_apply_and_returns_owned_output() -> None:
    operator = DiagonalOperator()
    state = np.array([4.0, 5.0])
    applied = operator.apply(state)
    multiplied = operator @ state

    np.testing.assert_array_equal(applied, np.array([8.0, -15.0]))
    np.testing.assert_array_equal(multiplied, applied)
    assert multiplied.flags.c_contiguous and multiplied.flags.owndata
    assert not np.shares_memory(multiplied, state)


def test_dense_inspection_is_derived_as_an_owned_c_contiguous_copy() -> None:
    operator = DiagonalOperator()
    dense = operator.to_dense()
    sparse = operator.matrix

    np.testing.assert_array_equal(dense, np.diag([2.0, -3.0]))
    assert dense.dtype == np.dtype(np.float64)
    assert dense.flags.c_contiguous and dense.flags.owndata
    dense[:] = 0.0
    sparse.data[:] = 0.0
    np.testing.assert_array_equal(operator.to_dense(), np.diag([2.0, -3.0]))


def test_real_scaling_preserves_endpoints_shape_dtype_and_observable_factor() -> None:
    operator = DiagonalOperator()
    scaled = operator.scaled(np.float32(-0.5))
    state = np.array([4.0, 5.0])
    result = scaled @ state

    assert scaled.domain is operator.domain
    assert scaled.codomain is operator.codomain
    assert scaled.shape == operator.shape
    assert scaled.dtype == np.dtype(np.float64)
    np.testing.assert_array_equal(result, np.array([-4.0, 7.5]))
    assert result.flags.c_contiguous and result.flags.owndata
    assert not np.shares_memory(result, state)


def test_complex_scaling_preserves_endpoints_and_canonical_complex_behavior() -> None:
    operator = DiagonalOperator()
    factor = np.complex64(1.0 + 2.0j)
    scaled = operator.scaled(factor)
    factor = np.complex64(99.0 + 99.0j)
    state = np.array([4.0, 5.0])
    result = scaled @ state

    assert scaled.domain is operator.domain
    assert scaled.codomain is operator.codomain
    assert scaled.shape == operator.shape
    assert scaled.dtype == np.dtype(np.complex128)
    np.testing.assert_array_equal(result, np.array([8.0 + 16.0j, -15.0 - 30.0j]))
    assert result.flags.c_contiguous and result.flags.owndata
    assert not np.shares_memory(result, state)
    assert factor == np.complex64(99.0 + 99.0j)


@pytest.mark.parametrize("factor", [True, np.bool_(False), "2", [2.0], object()])
def test_scaling_rejects_wrong_scalar_semantic_types(factor: object) -> None:
    with pytest.raises(TypeError):
        DiagonalOperator().scaled(factor)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "factor", [np.nan, np.inf, -np.inf, complex(np.nan, 0), complex(0, np.inf), 10**400]
)
def test_scaling_rejects_nonfinite_or_unrepresentable_numeric_factors(
    factor: object,
) -> None:
    with pytest.raises(ValueError):
        DiagonalOperator().scaled(factor)  # type: ignore[arg-type]


def test_scaled_operator_is_observably_immutable_and_has_no_public_composition() -> (
    None
):
    operator = DiagonalOperator()
    scaled = operator.scaled(2.0)

    for attribute, value in (("shape", (1, 1)), ("dtype", np.dtype(int)), ("extra", 1)):
        with pytest.raises(AttributeError):
            setattr(scaled, attribute, value)
    assert not hasattr(LinearOperator, "compose")
    assert not hasattr(DiscreteLinearOperator1D, "compose")
    assert not hasattr(operator, "compose")
    assert not hasattr(scaled, "compose")


def test_concrete_state_space_is_not_a_new_general_base_requirement() -> None:
    space = HomogeneousDirichletStateSpace1D(UniformGrid1D(0.0, 1.0, 5))

    assert not isinstance(space, LinearOperator)
    assert isinstance(DiagonalOperator().domain, tuple)
