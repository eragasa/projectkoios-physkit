# `DenseBlockHoppingBlochEvaluation` numerics

The result requires one dense complex matrix per coordinate row. Every matrix must
have shape `(basis_size, basis_size)` and the exact model energy unit. The coordinate
matrix must be unitless and have width equal to the lattice dimension. Empty
coordinate inventories correlate with empty matrix tuples.

The result performs structural correlation checks only; it does not recompute the
Fourier sum.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
