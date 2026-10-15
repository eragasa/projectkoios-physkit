# Class `BlockHoppingTerm`

## Current contract

`BlockHoppingTerm` binds one exact `LatticeDisplacement` to one finite immutable
nonempty square `ComplexMatrixQuantity`. Matrix rows and columns refer to the same
ordered primitive-cell basis identified by the owning `BlockHoppingModel`.

The record does not infer the negative-displacement partner or impose onsite
Hermiticity.

## Mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py::BlockHoppingTerm`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/test__BlockHoppingTerm__init.py`
