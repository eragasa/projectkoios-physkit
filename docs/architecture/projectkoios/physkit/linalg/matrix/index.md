# Package `projectkoios.physkit.linalg.matrix`

## Responsibility

This package groups operations whose primary represented object is a matrix. Each
mathematical operation occupies an explicit subpackage with its own representation,
algorithm-selection, result, and failure contracts.

## Public structure

- [`norm`](norm/index.md) owns quantity-aware matrix-norm evaluation.

The package `__init__.py` exports no identities. Callers import the norm integration
surface directly from `projectkoios.physkit.linalg.matrix.norm`. Neither `matrix` nor
`linalg` flattens norm classes into a broader namespace.

## Ownership boundary

`projectkoios.physkit.linalg.eigen` owns the established eigensystem and spectral-
decomposition contracts. The `matrix` package does not alias, wrap, or relocate those
contracts. No matrix-decomposition package is part of this package contract.

`projectkoios.physkit.linalg.matrices` remains the distinct module for matrix
verification and raw-to-quantity canonicalization. The singular `matrix` package owns
mathematical operation families; the plural `matrices` module owns representation
verification.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/matrix/__init__.py`
- Tests: `tests/projectkoios/physkit/linalg/matrix/test__matrix_public_api.py`
- Norm package: `src/python/projectkoios/physkit/linalg/matrix/norm/`

## Related pages

- [Package `projectkoios.physkit.linalg`](../index.md)
