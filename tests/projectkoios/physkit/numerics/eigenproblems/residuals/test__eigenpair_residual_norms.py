"""Verification of unit-free eigenpair residual norms."""

import warnings

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.numerics.eigenproblems.residuals import (
    eigenpair_residual_norms,
)


class TestEigenpairResidualNorms:
    """Verify generic dense and sparse numerical residual evaluation."""

    def test_returns_zero_for_exact_real_eigenpairs(self) -> None:
        operator = np.diag(np.array([2.0, 5.0], dtype=np.float64))
        eigenvalues = np.array([2.0, 5.0], dtype=np.float64)
        eigenvectors = np.eye(2, dtype=np.float64)

        norms = eigenpair_residual_norms(operator, eigenvalues, eigenvectors)

        np.testing.assert_array_equal(norms, np.zeros(2, dtype=np.float64))
        assert norms.dtype == np.dtype(np.float64)

    def test_returns_zero_for_exact_complex_eigenpairs(self) -> None:
        operator = np.diag(np.array([1.0, 2.0], dtype=np.complex128))
        eigenvalues = np.array([1.0 + 0.0j, 2.0 + 0.0j], dtype=np.complex128)
        eigenvectors = np.array(
            [[1.0j, 0.0j], [0.0j, (1.0 + 1.0j) / np.sqrt(2.0)]],
            dtype=np.complex128,
        )

        dense_norms = eigenpair_residual_norms(
            operator,
            eigenvalues,
            eigenvectors,
        )
        sparse_norms = eigenpair_residual_norms(
            sparse.csr_array(operator),
            eigenvalues,
            eigenvectors,
        )

        np.testing.assert_array_equal(dense_norms, np.zeros(2, dtype=np.float64))
        np.testing.assert_array_equal(sparse_norms, dense_norms)

    def test_returns_hand_derived_nonzero_residual(self) -> None:
        operator = np.diag(np.array([2.0, 4.0], dtype=np.float64))
        eigenvalues = np.array([3.0], dtype=np.float64)
        eigenvectors = np.array([[1.0], [2.0]], dtype=np.float64)

        norms = eigenpair_residual_norms(operator, eigenvalues, eigenvectors)

        np.testing.assert_allclose(norms, [np.sqrt(5.0)], rtol=0.0, atol=0.0)

    @pytest.mark.parametrize(
        "format_name",
        ["csr", "csc", "coo", "lil", "dok", "dia", "bsr"],
    )
    def test_dense_and_sparse_formats_are_equivalent(
        self,
        format_name: str,
    ) -> None:
        dense = np.array([[2.0, -1.0], [-1.0, 3.0]], dtype=np.float64)
        eigenvalues = np.array([1.0, 4.0], dtype=np.float64)
        eigenvectors = np.array([[1.0, 0.5], [0.0, 1.0]], dtype=np.float64)
        sparse_type = getattr(sparse, f"{format_name}_array")

        dense_norms = eigenpair_residual_norms(
            dense,
            eigenvalues,
            eigenvectors,
        )
        sparse_norms = eigenpair_residual_norms(
            sparse_type(dense),
            eigenvalues,
            eigenvectors,
        )

        np.testing.assert_array_equal(sparse_norms, dense_norms)

    def test_sparse_operator_is_not_densified(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        operator = sparse.csr_array(np.diag(np.array([2.0, 5.0], dtype=np.float64)))

        def reject_dense_conversion(*args: object, **kwargs: object) -> None:
            del args, kwargs
            raise AssertionError("sparse residual evaluation must not densify")

        monkeypatch.setattr(sparse.csr_array, "toarray", reject_dense_conversion)

        norms = eigenpair_residual_norms(
            operator,
            np.array([2.0, 5.0], dtype=np.float64),
            np.eye(2, dtype=np.float64),
        )

        np.testing.assert_array_equal(norms, np.zeros(2, dtype=np.float64))

    @pytest.mark.parametrize(
        "operator",
        [
            np.eye(2, dtype=np.bool_),
            np.eye(2, dtype=np.int64),
            np.eye(2, dtype=np.float32),
            np.array([["1.0", "0.0"], ["0.0", "1.0"]]),
            np.array([[1.0, 0.0], [0.0, 1.0]], dtype=object),
            np.eye(2, dtype=np.complex64),
            sparse.eye(2, dtype=np.bool_, format="csr"),
            sparse.eye(2, dtype=np.int64, format="csr"),
            sparse.eye(2, dtype=np.float32, format="csr"),
            sparse.eye(2, dtype=np.complex64, format="csr"),
        ],
    )
    def test_rejects_coercive_operator_dtypes(self, operator: object) -> None:
        with pytest.raises(TypeError, match="exactly float64 or complex128"):
            eigenpair_residual_norms(  # type: ignore[arg-type]
                operator,
                np.ones(1, dtype=np.float64),
                np.ones((2, 1), dtype=np.float64),
            )

    @pytest.mark.parametrize(
        ("eigenvalues", "eigenvectors", "name"),
        [
            (
                np.array([True], dtype=np.bool_),
                np.ones((2, 1), dtype=np.float64),
                "eigenvalues",
            ),
            (
                np.array([1], dtype=np.int64),
                np.ones((2, 1), dtype=np.float64),
                "eigenvalues",
            ),
            (
                np.array([1.0], dtype=np.float32),
                np.ones((2, 1), dtype=np.float64),
                "eigenvalues",
            ),
            (
                np.array(["1.0"]),
                np.ones((2, 1), dtype=np.float64),
                "eigenvalues",
            ),
            (
                np.array([1.0], dtype=object),
                np.ones((2, 1), dtype=np.float64),
                "eigenvalues",
            ),
            (
                np.ones(1, dtype=np.float64),
                np.ones((2, 1), dtype=np.bool_),
                "eigenvectors",
            ),
            (
                np.ones(1, dtype=np.float64),
                np.ones((2, 1), dtype=np.int64),
                "eigenvectors",
            ),
            (
                np.ones(1, dtype=np.float64),
                np.ones((2, 1), dtype=np.float32),
                "eigenvectors",
            ),
            (
                np.ones(1, dtype=np.float64),
                np.array([["1.0"], ["0.0"]]),
                "eigenvectors",
            ),
            (
                np.ones(1, dtype=np.float64),
                np.ones((2, 1), dtype=object),
                "eigenvectors",
            ),
        ],
    )
    def test_rejects_coercive_eigenpair_dtypes(
        self,
        eigenvalues: np.ndarray,
        eigenvectors: np.ndarray,
        name: str,
    ) -> None:
        with pytest.raises(
            TypeError,
            match=rf"{name} dtype must be exactly float64 or complex128",
        ):
            eigenpair_residual_norms(
                np.eye(2, dtype=np.float64),
                eigenvalues,
                eigenvectors,
            )

    @pytest.mark.parametrize(
        ("operator", "eigenvalues", "eigenvectors", "message"),
        [
            (
                np.ones((2, 3), dtype=np.float64),
                np.ones(1, dtype=np.float64),
                np.ones((2, 1), dtype=np.float64),
                "nonempty and square",
            ),
            (
                np.eye(2, dtype=np.float64),
                np.ones((1, 1), dtype=np.float64),
                np.ones((2, 1), dtype=np.float64),
                "eigenvalues must be a vector",
            ),
            (
                np.eye(2, dtype=np.float64),
                np.ones(1, dtype=np.float64),
                np.ones(2, dtype=np.float64),
                "eigenvectors must be a matrix",
            ),
            (
                np.eye(2, dtype=np.float64),
                np.ones(2, dtype=np.float64),
                np.ones((2, 1), dtype=np.float64),
                "must match operator and eigenvalue dimensions",
            ),
            (
                np.eye(2, dtype=np.float64),
                np.ones(3, dtype=np.float64),
                np.ones((2, 3), dtype=np.float64),
                "between one and the operator dimension",
            ),
        ],
    )
    def test_rejects_inconsistent_shapes(
        self,
        operator: np.ndarray,
        eigenvalues: np.ndarray,
        eigenvectors: np.ndarray,
        message: str,
    ) -> None:
        with pytest.raises(ValueError, match=message):
            eigenpair_residual_norms(operator, eigenvalues, eigenvectors)

    @pytest.mark.parametrize(
        ("operator", "eigenvalues", "eigenvectors", "name"),
        [
            (
                np.array([[np.nan]], dtype=np.float64),
                np.ones(1, dtype=np.float64),
                np.ones((1, 1), dtype=np.float64),
                "operator",
            ),
            (
                sparse.csr_array([[np.inf]], dtype=np.float64),
                np.ones(1, dtype=np.float64),
                np.ones((1, 1), dtype=np.float64),
                "operator",
            ),
            (
                np.eye(1, dtype=np.float64),
                np.array([np.nan], dtype=np.float64),
                np.ones((1, 1), dtype=np.float64),
                "eigenvalues",
            ),
            (
                np.eye(1, dtype=np.float64),
                np.ones(1, dtype=np.float64),
                np.array([[np.inf + 0.0j]], dtype=np.complex128),
                "eigenvectors",
            ),
        ],
    )
    def test_rejects_nonfinite_inputs(
        self,
        operator: object,
        eigenvalues: np.ndarray,
        eigenvectors: np.ndarray,
        name: str,
    ) -> None:
        with pytest.raises(ValueError, match=rf"{name} must contain only finite"):
            eigenpair_residual_norms(  # type: ignore[arg-type]
                operator,
                eigenvalues,
                eigenvectors,
            )

    def test_reports_nonfinite_arithmetic_without_warning_failure(self) -> None:
        maximum = np.finfo(np.float64).max
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with pytest.raises(ValueError, match="calculation must remain finite"):
                eigenpair_residual_norms(
                    np.array([[maximum]], dtype=np.float64),
                    np.zeros(1, dtype=np.float64),
                    np.array([[2.0]], dtype=np.float64),
                )
