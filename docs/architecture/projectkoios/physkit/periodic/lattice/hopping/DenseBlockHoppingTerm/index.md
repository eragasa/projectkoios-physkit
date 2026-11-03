# Class `DenseBlockHoppingTerm`

## Current contract

`DenseBlockHoppingTerm` binds one exact `LatticeDisplacement` to one finite immutable
nonempty square `ComplexMatrixQuantity`. Matrix rows and columns refer to the same
ordered primitive-cell basis identified by the owning `DenseBlockHoppingModel`.

The record does not infer the negative-displacement partner or impose onsite
Hermiticity.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py::DenseBlockHoppingTerm`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/test__DenseBlockHoppingTerm__init.py`
