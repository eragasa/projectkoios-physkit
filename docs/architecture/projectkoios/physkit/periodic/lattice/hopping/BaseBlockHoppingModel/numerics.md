# `BaseBlockHoppingModel` numerics

The abstract owner validates nonempty string identities, exact lattice dimension,
nonempty tuple ownership, homogeneous exact term type, common coefficient shape,
common unit, matching term dimensions, and lexicographically sorted unique
displacements. `basis_size` is derived from the first validated term.

Concrete classes supply their owned dense or sparse term type. The base does not
perform Fourier evaluation or matrix arithmetic.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
