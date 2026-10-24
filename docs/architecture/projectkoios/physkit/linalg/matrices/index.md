# Module `projectkoios.physkit.linalg.matrices`

## Current responsibility

This module owns nominal unions and actions for finite real or complex matrices
represented as raw NumPy dense arrays, raw SciPy CSR matrices, or canonical
PhysKit dense/CSR quantities. Quantity inputs are retained by identity. Raw
binary64 or complex128 inputs are copied into unitless PhysKit quantities.

`MatrixVerifier.is_hermitian` delegates to the unit-free sparse-safe numerical
predicate and never densifies sparse input. It establishes a represented
software property; it does not retain scientific acceptance policy.
Hermiticity is an operator-only property and is distinct from eigenpair residual
analysis of a correlated operator, eigenvalue vector, and eigenvector matrix,
which belongs to `linalg.eigen`.

## Public contracts

- `MatrixInput` is the supported raw-or-quantity input union.
- `VerifiedMatrix` is the canonical quantity union.
- [`MatrixVerifier`](MatrixVerifier/index.md) owns canonical input retention and
  matrix-property queries.
- [`VerifiedMatrixUtils`](VerifiedMatrixUtils/index.md) owns uniform explicit
  dense conversion.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/matrices.py`
- Tests: `tests/projectkoios/physkit/linalg/matrices/`
- Numerical kernel:
  `src/python/projectkoios/physkit/numerics/linear_algebra/matrix_properties.py`

## Navigation

- [Parent package](../index.md)
