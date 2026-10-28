# Internal dense complex-Hermitian spectral implementation

## Current responsibility

`_complex_hermitian.py` owns the dense complex
`SpectralDecompositionABC` product. Its constructor retains the exact
`ComplexMatrixQuantity`, validates a nonempty square exactly Hermitian
representation, passes only its unchanged numerical magnitude to the unit-free
dense eigensolver, computes complex eigenpair residuals, restores the original
matrix unit, and retains one maximum-rank public `Eigensystem`.

Construction is eager and computes all N eigenpairs. Eigenvalues and residual
norms retain the matrix unit. Eigenvectors are stored in a unitless
`ComplexMatrixQuantity`. The numerical backend canonicalizes each eigenvector's
global phase by making its largest-magnitude component real and nonnegative.
This phase convention does not select a unique basis within a degenerate
eigenspace.

After construction, `eigensystem(...)` is selection-only. `None`, zero, and N
return the retained full object by identity. A positive rank below N copies the
lowest eigenvalue, eigenvector, and residual prefixes into a new immutable
`Eigensystem`, retains the exact matrix object, and derives the maximum residual
from stored evidence. It neither solves again nor recomputes residuals.

This implementation deliberately accepts dense complex matrix quantities only.
Complex sparse matrices remain unsupported at the factory boundary rather than
being silently materialized as dense arrays. No tolerance, pass/fail policy, or
scientific validation claim is stored.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/eigen/_complex_hermitian.py`
- Factory tests: `tests/projectkoios/physkit/linalg/eigen/test__SpectralDecompositionFactory__create.py`
- Evidence and selection tests: `tests/projectkoios/physkit/linalg/eigen/test__ComplexHermitianSpectralDecomposition__eigensystem.py`
- Record tests: `tests/projectkoios/physkit/linalg/eigen/test__Eigensystem__init__.py`

## Navigation

- [Package](../index.md)
