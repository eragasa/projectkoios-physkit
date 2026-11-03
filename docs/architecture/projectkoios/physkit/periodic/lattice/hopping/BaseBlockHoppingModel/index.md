# Class `BaseBlockHoppingModel`

## Current contract

`BaseBlockHoppingModel` is the abstract behavioral owner shared by dense and sparse
block-hopping models. Concrete models provide identity, dimension, homogeneous term
inventory, energy unit, energy reference, basis identity, and exact storage kind.
The base derives `basis_size` and enforces common identity, dimension, shape, unit,
displacement-order, and uniqueness invariants. Concrete models supply the exact
owned term class, so sparse terms cannot enter a dense model and dense terms cannot
enter a sparse model.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py::BaseBlockHoppingModel`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/test__BaseBlockHoppingModel__inheritance.py`
