"""Verification of the unit-free Hermiticity predicate."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.numerics.linear_algebra.matrix_properties import (
    matrix_is_hermitian,
)


class TestMatrixIsHermitian:
    """Verify dense and sparse Hermiticity without model policy."""

    def test_accepts_real_symmetric_and_complex_hermitian_matrices(self) -> None:
        real = np.array([[2.0, -1.0], [-1.0, 3.0]], dtype=np.float64)
        complex_matrix = np.array(
            [[2.0 + 0.0j, 1.0 + 2.0j], [1.0 - 2.0j, 3.0 + 0.0j]],
            dtype=np.complex128,
        )

        assert matrix_is_hermitian(real)
        assert matrix_is_hermitian(complex_matrix)
        assert not matrix_is_hermitian(real[:, :1])

    def test_checks_sparse_input_without_densifying(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        matrix = sparse.csr_array(
            np.array(
                [[2.0 + 0.0j, 1.0 + 2.0j], [1.0 - 2.0j, 3.0 + 0.0j]],
                dtype=np.complex128,
            )
        )

        def reject_dense_conversion(*args: object, **kwargs: object) -> None:
            del args, kwargs
            raise AssertionError("sparse verification must not densify")

        monkeypatch.setattr(sparse.csr_array, "toarray", reject_dense_conversion)

        assert matrix_is_hermitian(matrix)

    def test_applies_explicit_tolerances(self) -> None:
        matrix = np.array([[1.0, 1.0 + 1.0e-10], [1.0, 2.0]], dtype=np.float64)

        assert not matrix_is_hermitian(matrix)
        assert matrix_is_hermitian(matrix, atol=2.0e-10)

    @pytest.mark.parametrize("tolerance", [True, 0, np.float64(0.0)])
    def test_rejects_non_builtin_float_tolerances(self, tolerance: object) -> None:
        with pytest.raises(TypeError, match="built-in float"):
            matrix_is_hermitian(
                np.eye(2, dtype=np.float64),
                atol=tolerance,  # type: ignore[arg-type]
            )

    def test_rejects_nonfinite_or_unsupported_matrices(self) -> None:
        with pytest.raises(TypeError, match="float64 or complex128"):
            matrix_is_hermitian(np.eye(2, dtype=np.int64))
        with pytest.raises(ValueError, match="finite"):
            matrix_is_hermitian(
                np.array([[1.0, np.inf], [np.inf, 1.0]], dtype=np.float64)
            )
