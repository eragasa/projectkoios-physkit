# Internal real-symmetric spectral implementation

## Current responsibility

`_real_symmetric.py` owns the current concrete `SpectralDecompositionABC`
product. Its constructor retains the exact real dense/CSR matrix quantity,
validates exact symmetry without densifying, passes only a raw magnitude/CSR
representation to the unit-free complete eigensolver, computes residuals with
`eigenpair_residual_norms`, restores units, and retains one maximum-rank public
`Eigensystem`.

Construction is deliberately eager and may be expensive. Complete sparse
matrices retain the backend's exact-tridiagonal restriction. The represented
implementation never calls the backend's `execute_lowest`; that raw numerical
capability remains available to other numerical users.

After construction, `eigensystem(...)` is selection-only. `None`, zero, and N
return the retained full object by identity. A positive rank below N copies the
lowest eigenvalue/eigenvector/residual prefixes into a new immutable
`Eigensystem`, retains the exact matrix object, and derives the maximum residual
from the stored residual prefix. It neither solves again nor recomputes
residuals.

The eigensystem uses the matrix unit for eigenvalues and residual evidence and
Unitless real eigenvectors. No tolerance or pass/fail policy is stored.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/eigen/_real_symmetric.py`
- Construction tests: `tests/projectkoios/physkit/linalg/eigen/test__SpectralDecomposition__init__.py`
- Selection tests: `tests/projectkoios/physkit/linalg/eigen/test__SpectralDecomposition__eigensystem.py`
- Record tests: `tests/projectkoios/physkit/linalg/eigen/test__Eigensystem__init__.py`

## Navigation

- [Package](../index.md)
