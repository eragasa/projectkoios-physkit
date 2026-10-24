# Module `projectkoios.physkit.linalg.eigen.factory`

## Current responsibility

This module owns the top-level, non-nested `SpectralDecompositionFactory` used
by the public façade. The factory performs bounded nominal dispatch and concrete
construction. Concrete construction is intentionally eager: selecting the real
product triggers its complete N-state eigensolve and residual computation before
`create(...)` returns. The factory does not own the numerical algorithm or any
tolerance policy.

Real dense/CSR quantities select the internal real-symmetric implementation.
The eager complete sparse path retains the exact-tridiagonal numerical backend
restriction. Malformed real or complex matrices reject. Exact
complex-Hermitian quantities raise `UnsupportedSpectralDecompositionError`
because no complex implementation exists.

The factory is a package collaboration seam, not a supported package export.
No plugin registry or dynamic algorithm registration is provided.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/eigen/factory.py`
- Tests: `tests/projectkoios/physkit/linalg/eigen/test__SpectralDecompositionFactory__create.py`

## Navigation

- [Package](../index.md)
