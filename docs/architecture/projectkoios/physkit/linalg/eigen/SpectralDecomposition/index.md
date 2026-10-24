# `SpectralDecomposition`

## Responsibility

`SpectralDecomposition(matrix)` is the sole supported represented construction
façade. It invokes the separate top-level `SpectralDecompositionFactory` to
select and eagerly construct a `SpectralDecompositionABC` product. Construction
performs the complete N-state eigensolve and residual computation before the
façade returns. This intentionally permits expensive hidden constructor work.
The complete sparse real path retains the numerical backend's exact
tridiagonal restriction, so unsupported nontridiagonal sparse matrices fail at
construction.

`matrix` returns the exact input quantity retained by the selected product.
`max_rank` returns the matrix dimension N, not a numerical rank estimate.
`eigensystem(n_rank=None)` performs no numerical decomposition:

- `None`, zero, or N returns the retained maximum-rank `Eigensystem` by identity;
- a positive built-in integer below N returns an immutable lowest-prefix
  `Eigensystem` whose matrix is the exact original object;
- booleans, NumPy integers, floats, negative values, and values above N reject.

Raw NumPy/SciPy values are not accepted. No tolerance, relative scaling, or
pass/fail policy belongs to decomposition or retained selection.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/eigen/facade.py::SpectralDecomposition`
- Construction tests: `tests/projectkoios/physkit/linalg/eigen/test__SpectralDecomposition__init__.py`
- Selection tests: `tests/projectkoios/physkit/linalg/eigen/test__SpectralDecomposition__eigensystem.py`

## Navigation

- [Package](../index.md)
