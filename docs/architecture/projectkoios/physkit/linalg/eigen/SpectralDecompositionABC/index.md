# `SpectralDecompositionABC`

## Responsibility

`SpectralDecompositionABC` is the nominal product contract for every concrete
represented decomposition. A concrete product must finish its complete
maximum-rank eigensolve and residual computation before construction returns.
It exposes the exact retained matrix, dimension-only `max_rank`, and
`eigensystem(n_rank=None)` selection over already-retained evidence.

The selection method must not invoke an eigensolver or recompute residuals.
`None`, zero, and `max_rank` denote the same retained full eigensystem. A
positive smaller rank denotes the lowest retained prefix.

The ABC does not dispatch matrix families, coerce representations, choose a
numerical algorithm, define tolerance policy, or estimate numerical matrix
rank.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/eigen/base.py::SpectralDecompositionABC`
- Factory tests: `tests/projectkoios/physkit/linalg/eigen/test__SpectralDecompositionFactory__create.py`
- Selection tests: `tests/projectkoios/physkit/linalg/eigen/test__SpectralDecomposition__eigensystem.py`

## Navigation

- [Package](../index.md)
