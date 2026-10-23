"""Unit-free numerical predicates for finite dense and sparse matrices."""

from __future__ import annotations

import math

import numpy as np
from scipy import linalg, sparse

from projectkoios.physkit.numerics.typing.numpy.arrays import (
    ComplexMatrix,
    RealMatrix,
)
from projectkoios.physkit.numerics.typing.scipy.sparse import SparseMatrix

type HermiticityMatrix = RealMatrix | ComplexMatrix | SparseMatrix


def matrix_is_hermitian(
    matrix: HermiticityMatrix,
    *,
    rtol: float = 0.0,
    atol: float = 0.0,
) -> bool:
    """Return whether one finite real or complex matrix is Hermitian.

    Dense inputs use SciPy's Hermiticity predicate. Sparse inputs compare the
    represented matrix with its conjugate transpose without materializing a
    dense matrix. Tolerances follow NumPy ``allclose`` elementwise semantics.
    """
    _require_tolerance(rtol, "rtol")
    _require_tolerance(atol, "atol")
    if isinstance(matrix, np.ndarray):
        _require_dtype(matrix.dtype)
        if matrix.ndim != 2:
            raise ValueError("matrix must be two-dimensional")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("matrix must contain only finite values")
        rows, columns = matrix.shape
        if rows != columns:
            return False
        return bool(linalg.ishermitian(matrix, rtol=rtol, atol=atol))
    if not sparse.issparse(matrix):
        raise TypeError("matrix must be a NumPy or SciPy sparse matrix")
    _require_dtype(matrix.dtype)
    represented = sparse.csr_array(matrix)
    represented.sum_duplicates()
    represented.eliminate_zeros()
    represented.sort_indices()
    if not np.all(np.isfinite(represented.data)):
        raise ValueError("matrix must contain only finite values")
    rows, columns = represented.shape
    if rows != columns:
        return False
    if rows == 0:
        return True
    adjoint = represented.conjugate().T.tocsr()
    difference = sparse.coo_array(represented - adjoint)
    difference.sum_duplicates()
    difference.eliminate_zeros()
    if difference.nnz == 0:
        return True
    reference = np.asarray(adjoint[difference.row, difference.col]).reshape(-1)
    return bool(np.all(np.abs(difference.data) <= atol + rtol * np.abs(reference)))


def _require_dtype(dtype: np.dtype[np.generic]) -> None:
    """Require exact binary64 real or complex storage."""
    if dtype not in (np.dtype(np.float64), np.dtype(np.complex128)):
        raise TypeError("matrix dtype must be exactly float64 or complex128")


def _require_tolerance(value: float, name: str) -> None:
    """Require one finite nonnegative built-in floating tolerance."""
    if type(value) is not float:
        raise TypeError(f"{name} must be a built-in float")
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
