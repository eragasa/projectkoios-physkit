# Module `projectkoios.physkit.linalg.matrix_norms.base`

## Responsibility

The module owns the common quantity-aware matrix-norm contract and the public semantic
dispatch facade.

## Owned classes

- [`BaseMatrixQuantityNorm`](BaseMatrixQuantityNorm/index.md) owns common quantity
  retention, representation extraction, and exact unit restoration.
- [`MatrixNorm`](MatrixNorm/index.md) selects one concrete quantity-aware norm from the
  exact canonical `MatrixNormEnum`.

The module imports Optimization's exact `MatrixNormABC` for implementation typing and
its exact `MatrixNormEnum` for public semantics. It does not redefine or re-export the
unit-free ABC through the PhysKit package facade.

## Exclusions

The module owns no numerical kernel, raw numerical facade, string selector, factory,
static utility behavior, fallback, result cache, or scientific acceptance policy.

## Related pages

- [Matrix-norm package](../index.md)
- [Package `projectkoios.physkit.linalg`](../../index.md)
