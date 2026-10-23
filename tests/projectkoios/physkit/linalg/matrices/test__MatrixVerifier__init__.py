"""Verification of ``MatrixVerifier`` construction."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.linalg import MatrixVerifier
from projectkoios.physkit.units import (
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    SparseMatrixQuantity,
    Unitless,
)


class TestMatrixVerifierInit:
    """Verify canonical immutable matrix retention."""

    def test_retains_quantity_input_by_identity(self) -> None:
        matrix = ComplexSparseMatrixQuantity.from_csr(
            sparse.eye(2, dtype=np.complex128, format="csr"),
            Unitless(),
        )

        verifier = MatrixVerifier(matrix)

        assert verifier.matrix is matrix

    def test_copies_raw_dense_input_into_immutable_quantity(self) -> None:
        source = np.eye(2, dtype=np.float64)

        verifier = MatrixVerifier(source)
        source[:] = 0.0

        assert isinstance(verifier.matrix, MatrixQuantity)
        np.testing.assert_array_equal(verifier.matrix.magnitude, np.eye(2))
        assert not verifier.matrix.magnitude.flags.writeable

    def test_copies_raw_csr_input_into_immutable_quantity(self) -> None:
        source = sparse.csr_array(np.eye(2, dtype=np.float64))

        verifier = MatrixVerifier(source)
        source.data[:] = 0.0

        assert isinstance(verifier.matrix, SparseMatrixQuantity)
        np.testing.assert_array_equal(
            verifier.matrix.to_dense().magnitude,
            np.eye(2),
        )

    def test_rejects_unsupported_raw_dtype_and_sparse_format(self) -> None:
        with pytest.raises(TypeError, match="float64 or complex128"):
            MatrixVerifier(np.eye(2, dtype=np.int64))
        with pytest.raises(TypeError, match="CSR-sparse"):
            MatrixVerifier(sparse.coo_array(np.eye(2, dtype=np.float64)))  # type: ignore[arg-type]
