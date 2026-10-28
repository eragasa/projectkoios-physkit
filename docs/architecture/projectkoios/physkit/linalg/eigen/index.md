# Package `projectkoios.physkit.linalg.eigen`

## Current responsibility

This package is the supported unit-aware spectral-decomposition boundary. The
public `SpectralDecomposition(matrix)` façade accepts only canonical PhysKit
dense/CSR matrix quantities and eagerly computes the complete N-state
spectral decomposition and residual evidence during construction. This eager
work is intentionally potentially expensive.

After construction, `eigensystem(n_rank=None)` only selects retained evidence.
`None`, zero, and N return the retained full eigensystem by identity. A positive
built-in integer below N returns an immutable lowest-prefix eigensystem without
calling an eigensolver or recomputing residuals. `max_rank` means exactly matrix
dimension N, never a tolerance-dependent numerical matrix rank.

Complete real-symmetric dense/CSR decomposition and complete dense
complex-Hermitian decomposition are implemented. Complex sparse matrices are
recognized and checked for exact Hermiticity but remain explicitly unsupported;
they are never silently densified.

## Public contracts

- [`SpectralDecomposition`](SpectralDecomposition/index.md) is the sole supported
  matrix-to-eager-decomposition façade.
- [`SpectralDecompositionABC`](SpectralDecompositionABC/index.md) is the nominal
  product contract implemented by every eager concrete decomposition.
- [`Eigensystem`](Eigensystem/index.md) is the public immutable evidence record.
- `UnsupportedSpectralDecompositionError` distinguishes a valid but unimplemented
  matrix family from malformed input.

`SpectralDecompositionFactory` is a top-level collaborator in `factory.py`, not
a nested class, but is intentionally omitted from package exports. The concrete
real decomposition identity is also an implementation detail.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/eigen/`
- Tests: `tests/projectkoios/physkit/linalg/eigen/`
- Unit-free eigensolvers:
  `projectkoios.physkit.numerics.eigenproblems.real_symmetric` and
  `projectkoios.physkit.numerics.eigenproblems.complex_hermitian`
- Unit-free residual kernel: `projectkoios.physkit.numerics.eigenproblems.residuals`

## Navigation

- [Factory module](factory/index.md)
- [Real implementation](real_symmetric/index.md)
- [Dense complex implementation](complex_hermitian/index.md)
- [Parent package](../index.md)
