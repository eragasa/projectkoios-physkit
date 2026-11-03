# `BlockHoppingModel` numerics

`create(...)` requires a nonempty tuple containing one exact concrete term type.
A homogeneous dense inventory selects `DenseBlockHoppingModel`; a homogeneous
sparse inventory selects `SparseBlockHoppingModel`. Mixed inventories and unsupported
values are rejected before concrete-model construction.

The selected concrete model performs the complete shared invariant validation. The
factory does not convert coefficient storage or reorder terms.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
