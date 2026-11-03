# `BaseBlockHoppingTerm` numerics

The abstract owner requires an exact `LatticeDisplacement`, a two-dimensional
coefficient shape, an explicit unit, and an exact storage kind. Shared validation
requires a nonzero row count and equal row and column counts. `basis_size` is derived
as the first coefficient-shape component.

The base performs no coefficient arithmetic and cannot be constructed directly.

## Numerical validation behavior

None. The abstract term owner performs no numerical arithmetic.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
