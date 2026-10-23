"""Verification of ``MatrixVerifier.is_hermitian``."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.linalg import MatrixVerifier
from projectkoios.physkit.units import ComplexSparseMatrixQuantity, Unitless


class TestMatrixVerifierIsHermitian:
    """Verify quantity-aware dense and sparse decisions."""

    def test_checks_sparse_quantity_without_densifying(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        matrix = ComplexSparseMatrixQuantity.from_csr(
            sparse.csr_array(
                np.array(
                    [[1.0 + 0.0j, 2.0j], [-2.0j, 3.0 + 0.0j]],
                    dtype=np.complex128,
                )
            ),
            Unitless(),
        )

        def reject_dense_conversion(*args: object, **kwargs: object) -> None:
            del args, kwargs
            raise AssertionError("sparse verification must not densify")

        monkeypatch.setattr(sparse.csr_array, "toarray", reject_dense_conversion)

        assert MatrixVerifier(matrix).is_hermitian()

    def test_forwards_strict_tolerance_validation(self) -> None:
        verifier = MatrixVerifier(np.eye(2, dtype=np.float64))

        with pytest.raises(TypeError, match="built-in float"):
            verifier.is_hermitian(rtol=True)  # type: ignore[arg-type]
        with pytest.raises(ValueError, match="finite and nonnegative"):
            verifier.is_hermitian(atol=-1.0)
