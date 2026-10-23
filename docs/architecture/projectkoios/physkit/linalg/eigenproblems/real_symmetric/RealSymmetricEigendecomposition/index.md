# `RealSymmetricEigendecomposition`

## Responsibility

`RealSymmetricEigendecomposition` is an immutable `ResultsObject` correlating
one exact real-symmetric `MatrixQuantity` or `SparseMatrixQuantity` with ordered
eigenvalues in the operator unit and unitless orthonormal eigenvector columns.
It may represent either a complete spectrum or an ordered lowest subset.

The result exposes `full_dimension`, `selected_count`, and `is_complete`. It
does not retain convergence policy, residual acceptance, model reduction, or
complex-Hermitian claims.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/eigenproblems/real_symmetric.py::RealSymmetricEigendecomposition`
- Tests: `tests/projectkoios/physkit/linalg/eigenproblems/real_symmetric/test__RealSymmetricEigendecomposition__init__.py`

## Navigation

- [Module](../index.md)
