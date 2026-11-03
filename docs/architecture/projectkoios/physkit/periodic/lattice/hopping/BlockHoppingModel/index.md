# Class `BlockHoppingModel`

## Current contract

`BlockHoppingModel` is the instantiable construction action for block-hopping models.
Its instance `create(...)` method selects `DenseBlockHoppingModel` for a homogeneous dense-term
inventory and `SparseBlockHoppingModel` for a homogeneous sparse-term inventory.
Mixed and empty inventories are rejected. Both products retain the common
`BaseBlockHoppingModel` abstraction.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py::BlockHoppingModel`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/test__BlockHoppingFactories__create.py`
