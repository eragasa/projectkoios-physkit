# Package `projectkoios.physkit.linalg`

## Current responsibility

This package is the public generic linear-algebra boundary for finite raw and
unit-aware matrix representations. It composes canonical quantities from
`projectkoios.physkit.units` with unit-free algorithms from
`projectkoios.physkit.numerics` without duplicating either ownership layer.

Callers may verify raw NumPy dense or SciPy CSR matrices and their PhysKit
quantity equivalents. Quantity-aware matrix norms accept only canonical real or complex
dense or sparse quantities, delegate raw magnitudes to exact Optimization owners, and
restore the exact input unit. The explicit `linalg.eigen` namespace accepts unit-aware
dense or sparse matrices, restores units around unit-free backends, and retains
per-eigenpair residual evidence without acceptance policy.

## Public contract

The package facade re-exports deliberate exact identities from `matrices` and
`matrix_norms`. Matrix-norm semantics and iterative controls reuse their exact
Optimization identities; PhysKit adds only quantity ownership, explicit sparse spectral
selection, and unit restoration. Spectral-decomposition contracts remain under the
explicit `linalg.eigen` namespace. Only real-symmetric spectral decomposition is
currently established. Generic matrix verification and matrix norms support real and
complex matrices but do not imply a validated complex-Hermitian decomposition contract.

## Navigation

- [Module `matrices`](matrices/index.md)
- [Package `matrix_norms`](matrix_norms/index.md)
- [Package `eigen`](eigen/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/`
- Tests: `tests/projectkoios/physkit/linalg/`
