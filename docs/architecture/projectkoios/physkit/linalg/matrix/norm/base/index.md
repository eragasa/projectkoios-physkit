# Module `projectkoios.physkit.linalg.matrix.norm.base`

## Responsibility

This module owns the shared quantity-aware matrix-norm contract and the public semantic
dispatch facade. It contains no numerical norm kernel.

## Owned classes

- [`BaseMatrixQuantityNorm`](BaseMatrixQuantityNorm/index.md) owns nominal quantity
  retention, dense/CSR representation access, and exact unit restoration.
- [`MatrixNorm`](MatrixNorm/index.md) maps one exact `MatrixNormEnum` member to one
  concrete quantity-aware evaluator.

## External identities

The module imports Optimization's exact `MatrixNormEnum`. It neither redefines nor wraps
that enum. Optimization's `MatrixNormABC` remains the separate unit-free result contract
and is not re-exported by PhysKit.

## Internal state

`MatrixNorm` installs two immutable `MappingProxyType` mappings per facade instance:

- semantic enum to stable diagnostic name;
- semantic enum to concrete quantity-aware evaluator class.

Both mappings use underscore names, have no public accessor, and are excluded from the
facade representation and equality semantics. They support deterministic dispatch only.

## Exclusions

The module performs no numerical reduction, singular-value decomposition, unit
conversion, sparse materialization, retry, fallback, result caching, thresholding, or
scientific acceptance.

## Related pages

- [Matrix-norm package](../index.md)
- [Package `projectkoios.physkit.linalg.matrix`](../../index.md)
- [Package `projectkoios.physkit.linalg`](../../../index.md)
