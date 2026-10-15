"""Verification of complete unit-free real-symmetric eigensolution."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.numerics.eigenproblems.real_symmetric import (
    RealSymmetricEigenpairSolver,
)
from projectkoios.physkit.units import MatrixQuantity, Unitless


def test_returns_ordered_reconstructing_dense_eigenpairs() -> None:
    operator = np.array([[2.0, -1.0], [-1.0, 2.0]])
    result = RealSymmetricEigenpairSolver().execute(operator)

    np.testing.assert_allclose(result.eigenvalues, [1.0, 3.0])
    np.testing.assert_allclose(
        result.eigenvectors @ np.diag(result.eigenvalues) @ result.eigenvectors.T,
        operator,
        atol=8.0 * np.finfo(np.float64).eps,
    )


def test_solves_sparse_tridiagonal_without_densifying(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    operator = sparse.csr_array(
        sparse.diags(
            ([-1.0, -1.0], [2.0, 2.0, 2.0], [-1.0, -1.0]),
            offsets=(-1, 0, 1),
            format="csr",
        )
    )

    def reject_dense_conversion(*args: object, **kwargs: object) -> None:
        del args, kwargs
        raise AssertionError("sparse input must not be materialized as dense")

    monkeypatch.setattr(sparse.csr_array, "toarray", reject_dense_conversion)
    result = RealSymmetricEigenpairSolver().execute(operator)

    reconstructed = (
        result.eigenvectors @ np.diag(result.eigenvalues) @ result.eigenvectors.T
    )
    expected = np.array([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]])
    np.testing.assert_allclose(
        reconstructed,
        expected,
        atol=16.0 * np.finfo(np.float64).eps,
    )


def test_rejects_unit_bearing_operator_at_numerical_boundary() -> None:
    with pytest.raises(TypeError, match="NumPy or SciPy"):
        RealSymmetricEigenpairSolver().execute(
            MatrixQuantity(np.eye(2), Unitless())  # type: ignore[arg-type]
        )


def test_rejects_nonsymmetric_and_nontridiagonal_sparse_operators() -> None:
    solver = RealSymmetricEigenpairSolver()
    with pytest.raises(ValueError, match="exactly real symmetric"):
        solver.execute(np.array([[1.0, 1.0], [0.0, 1.0]]))
    nontridiagonal = sparse.csr_array(
        [[2.0, -1.0, 0.5], [-1.0, 2.0, -1.0], [0.5, -1.0, 2.0]]
    )
    with pytest.raises(ValueError, match="tridiagonal"):
        solver.execute(nontridiagonal)
