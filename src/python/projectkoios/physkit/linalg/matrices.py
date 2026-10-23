"""Quantity-aware finite dense- and sparse-matrix verification."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from projectkoios.physkit.numerics.linear_algebra.matrix_properties import (
    matrix_is_hermitian,
)
from projectkoios.physkit.numerics.typing.numpy.arrays import (
    ComplexMatrix,
    RealMatrix,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    SparseMatrixQuantity,
    Unitless,
)

type CsrSparseMatrix = sparse.csr_matrix | sparse.csr_array
type VerifiedMatrix = (
    MatrixQuantity
    | ComplexMatrixQuantity
    | SparseMatrixQuantity
    | ComplexSparseMatrixQuantity
)
type MatrixInput = RealMatrix | ComplexMatrix | CsrSparseMatrix | VerifiedMatrix


class VerifiedMatrixUtils:
    """Apply representation conversions to one verified matrix quantity."""

    __slots__ = ("matrix",)

    matrix: VerifiedMatrix

    def __init__(self, matrix: VerifiedMatrix) -> None:
        """Retain one verified matrix representation."""
        if not isinstance(
            matrix,
            MatrixQuantity
            | ComplexMatrixQuantity
            | SparseMatrixQuantity
            | ComplexSparseMatrixQuantity,
        ):
            raise TypeError("matrix must be a verified matrix quantity")
        self.matrix = matrix

    def to_dense_matrix(self) -> MatrixQuantity | ComplexMatrixQuantity:
        """Return an immutable dense matrix quantity with the retained unit."""
        if isinstance(self.matrix, MatrixQuantity | ComplexMatrixQuantity):
            return self.matrix
        return self.matrix.to_dense()


class MatrixVerifier:
    """Verify properties of one immutable finite matrix representation.

    Quantity inputs are retained by identity. Raw dense arrays and CSR matrices
    must use exactly ``float64`` or ``complex128`` storage and are copied into a
    canonical unitless PhysKit quantity before verification.
    """

    __slots__ = ("matrix",)

    matrix: VerifiedMatrix

    def __init__(self, matrix: MatrixInput) -> None:
        """Validate and retain an immutable canonical representation."""
        if isinstance(
            matrix,
            MatrixQuantity
            | ComplexMatrixQuantity
            | SparseMatrixQuantity
            | ComplexSparseMatrixQuantity,
        ):
            self.matrix = matrix
            return
        if isinstance(matrix, np.ndarray):
            self._require_raw_dtype(matrix.dtype)
            if matrix.ndim != 2:
                raise ValueError("raw dense matrix must be two-dimensional")
            if not np.all(np.isfinite(matrix)):
                raise ValueError("raw dense matrix must contain only finite values")
            if matrix.dtype == np.dtype(np.complex128):
                self.matrix = ComplexMatrixQuantity(matrix, Unitless())
            else:
                self.matrix = MatrixQuantity(matrix, Unitless())
            return
        if isinstance(matrix, sparse.csr_matrix | sparse.csr_array):
            self._require_raw_dtype(matrix.dtype)
            if np.iscomplexobj(matrix.data):
                self.matrix = ComplexSparseMatrixQuantity.from_csr(
                    matrix,
                    Unitless(),
                )
            else:
                self.matrix = SparseMatrixQuantity.from_csr(matrix, Unitless())
            return
        raise TypeError(
            "matrix must be a dense or CSR-sparse real or complex matrix representation"
        )

    def is_hermitian(
        self,
        *,
        rtol: float = 0.0,
        atol: float = 0.0,
    ) -> bool:
        """Return whether the retained matrix is Hermitian without densifying CSR."""
        if isinstance(self.matrix, MatrixQuantity | ComplexMatrixQuantity):
            represented = self.matrix.magnitude
        else:
            represented = self.matrix.to_csr()
        return matrix_is_hermitian(represented, rtol=rtol, atol=atol)

    @staticmethod
    def _require_raw_dtype(dtype: np.dtype[np.generic]) -> None:
        """Require exact binary64 real or complex raw storage."""
        if dtype not in (np.dtype(np.float64), np.dtype(np.complex128)):
            raise TypeError("raw matrix dtype must be exactly float64 or complex128")
