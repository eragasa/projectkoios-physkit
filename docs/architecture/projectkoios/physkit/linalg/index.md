# Package `projectkoios.physkit.linalg`

## Current responsibility

This package is the public generic linear-algebra boundary for finite raw and
unit-aware matrix representations. It composes canonical quantities from
`projectkoios.physkit.units` with unit-free algorithms from
`projectkoios.physkit.numerics` without duplicating either ownership layer.

Callers may verify raw NumPy dense or SciPy CSR matrices and their PhysKit
quantity equivalents. Represented eigendecompositions accept unit-aware dense
or sparse matrices and restore units around the unit-free backend.

## Public contract

The package façade re-exports exact defining-module objects from `matrices` and
`eigenproblems`; it does not introduce wrapper subclasses or conversion copies.
Only real-symmetric eigendecomposition is currently established. Generic matrix
verification supports real and complex matrices but does not imply a validated
complex-Hermitian decomposition contract.

## Navigation

- [Module `matrices`](matrices/index.md)
- [Package `eigenproblems`](eigenproblems/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/`
- Tests: `tests/projectkoios/physkit/linalg/`
