# `SparseBlockHoppingBlochEvaluation` numerics

The result requires one canonical complex sparse matrix per coordinate row. Every
matrix must have shape `(basis_size, basis_size)` and the exact model energy unit.
The coordinate matrix must be unitless and match the lattice dimension. Empty
coordinate inventories correlate with empty matrix tuples.

The result does not densify matrices or recompute the Fourier sum.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
