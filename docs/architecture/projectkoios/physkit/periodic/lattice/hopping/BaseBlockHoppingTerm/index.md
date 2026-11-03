# Class `BaseBlockHoppingTerm`

## Current contract

`BaseBlockHoppingTerm` is the abstract behavioral owner shared by dense and sparse
matrix-valued hopping terms. Concrete terms must provide displacement,
`coefficient_shape`, energy unit, and exact storage kind. The base derives
`basis_size` and enforces the shared exact-displacement and nonempty-square-block
invariants. It defines no coefficient union and performs no representation
conversion.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py::BaseBlockHoppingTerm`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/test__BaseBlockHoppingModel__inheritance.py`
