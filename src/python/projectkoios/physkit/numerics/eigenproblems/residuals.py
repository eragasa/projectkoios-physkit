"""Unit-free numerical residuals for finite eigenpairs."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from projectkoios.physkit.numerics.typing.numpy.arrays import (
    ComplexMatrix,
    ComplexVector,
    RealMatrix,
    RealVector,
)
from projectkoios.physkit.numerics.typing.scipy.sparse import SparseMatrix

type EigenpairOperator = RealMatrix | ComplexMatrix | SparseMatrix
type EigenvalueArray = RealVector | ComplexVector
type EigenvectorArray = RealMatrix | ComplexMatrix

_SUPPORTED_DTYPES = (np.dtype(np.float64), np.dtype(np.complex128))


def eigenpair_residual_norms(
    operator: EigenpairOperator,
    eigenvalues: EigenvalueArray,
    eigenvectors: EigenvectorArray,
) -> RealVector:
    """Return one Euclidean residual norm per selected eigenvector column.

    The residual matrix is ``operator @ eigenvectors - eigenvectors *
    eigenvalues``. Inputs must already use exact ``float64`` or ``complex128``
    storage. Sparse operators remain sparse throughout multiplication.
    """
    represented_operator = _validate_operator(operator)
    rows, columns = represented_operator.shape
    if rows <= 0 or rows != columns:
        raise ValueError("operator must be nonempty and square")
    values = _validate_array(eigenvalues, "eigenvalues", dimensions=1)
    vectors = _validate_array(eigenvectors, "eigenvectors", dimensions=2)
    selected = values.size
    if selected <= 0 or selected > rows:
        raise ValueError(
            "eigenvalues must select between one and the operator dimension"
        )
    if vectors.shape != (rows, selected):
        raise ValueError("eigenvectors must match operator and eigenvalue dimensions")

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        residuals = represented_operator @ vectors - vectors * values[np.newaxis, :]
        if not np.all(np.isfinite(residuals)):
            raise ValueError("eigenpair residual calculation must remain finite")
        norms = np.linalg.norm(residuals, ord=2, axis=0)
    if not np.all(np.isfinite(norms)):
        raise ValueError("eigenpair residual norms must remain finite")
    return np.asarray(norms, dtype=np.float64)


def _validate_operator(
    operator: EigenpairOperator,
) -> RealMatrix | ComplexMatrix | sparse.csr_array:
    """Validate and return one dense or canonical-CSR numerical operator."""
    if isinstance(operator, np.ndarray):
        return _validate_array(operator, "operator", dimensions=2)
    if not sparse.issparse(operator):
        raise TypeError("operator must be a NumPy or SciPy sparse matrix")
    _validate_dtype(operator.dtype, "operator")
    if operator.ndim != 2:
        raise ValueError("operator must be a matrix")
    with np.errstate(over="ignore", invalid="ignore"):
        represented = sparse.csr_array(operator)
        represented.sum_duplicates()
        represented.eliminate_zeros()
        represented.sort_indices()
    if not np.all(np.isfinite(represented.data)):
        raise ValueError("operator must contain only finite values")
    return represented


def _validate_array(
    array: np.ndarray,
    name: str,
    *,
    dimensions: int,
) -> np.ndarray:
    """Validate one exact finite binary64 real or complex NumPy array."""
    if not isinstance(array, np.ndarray):
        raise TypeError(f"{name} must be a numpy.ndarray")
    _validate_dtype(array.dtype, name)
    if array.ndim != dimensions:
        expected = "a vector" if dimensions == 1 else "a matrix"
        raise ValueError(f"{name} must be {expected}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _validate_dtype(dtype: np.dtype[np.generic], name: str) -> None:
    """Require exact binary64 real or complex storage without coercion."""
    if dtype not in _SUPPORTED_DTYPES:
        raise TypeError(f"{name} dtype must be exactly float64 or complex128")
