"""Verification of ``VerifiedMatrixUtils.to_dense_matrix``."""

import numpy as np
from scipy import sparse

from projectkoios.physkit.linalg.matrices import VerifiedMatrixUtils
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    SparseMatrixQuantity,
)


class TestVerifiedMatrixUtilsToDenseMatrix:
    """Verify explicit dense materialization behavior."""

    def test_returns_dense_input_by_identity(self) -> None:
        matrix = MatrixQuantity(np.eye(2), PhysicalUnit("electron_volt"))

        assert VerifiedMatrixUtils(matrix).to_dense_matrix() is matrix

    def test_materializes_sparse_input_with_its_unit(self) -> None:
        unit = PhysicalUnit("electron_volt")
        matrix = SparseMatrixQuantity.from_csr(
            sparse.eye(2, dtype=np.float64, format="csr"),
            unit,
        )

        dense = VerifiedMatrixUtils(matrix).to_dense_matrix()

        assert isinstance(dense, MatrixQuantity)
        assert dense.unit == unit
        np.testing.assert_array_equal(dense.magnitude, np.eye(2))
