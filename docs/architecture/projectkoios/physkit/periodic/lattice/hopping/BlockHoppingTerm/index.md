# Class `BlockHoppingTerm`

## Current contract

`BlockHoppingTerm` is the instantiable construction action for block-hopping terms.
Its instance `create(...)` method selects `DenseBlockHoppingTerm` for an exact
`ComplexMatrixQuantity` and `SparseBlockHoppingTerm` for an exact
`ComplexSparseMatrixQuantity`. It performs no implicit dense–sparse conversion.
Both products retain the common `BaseBlockHoppingTerm` abstraction.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py::BlockHoppingTerm`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/test__BlockHoppingFactories__create.py`
