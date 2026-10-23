# Module `projectkoios.physkit.linalg.eigenproblems.real_symmetric`

## Current responsibility

This module owns represented real-symmetric eigendecompositions. The operator
and eigenvalues retain one exact unit, eigenvectors are unitless orthonormal
columns, and the numerical backend receives only stripped NumPy or SciPy
representations.

## Public classes

- [`RealSymmetricEigendecomposition`](RealSymmetricEigendecomposition/index.md)
- [`RealSymmetricEigendecomposer`](RealSymmetricEigendecomposer/index.md)

## Public type aliases

- `RealSymmetricOperator = MatrixQuantity | SparseMatrixQuantity`

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/eigenproblems/real_symmetric.py`
- Tests: `tests/projectkoios/physkit/linalg/eigenproblems/real_symmetric/`
- Numerical backend:
  `src/python/projectkoios/physkit/numerics/eigenproblems/real_symmetric.py`

## Navigation

- [Parent package](../index.md)
