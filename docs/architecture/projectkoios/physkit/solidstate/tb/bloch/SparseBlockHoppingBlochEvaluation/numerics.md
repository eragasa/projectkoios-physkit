# `SparseBlockHoppingBlochEvaluation` numerics

The result requires one canonical complex sparse matrix per coordinate row. Every
matrix must have shape `(basis_size, basis_size)` and the exact model energy unit.
The coordinate matrix must be unitless and match the lattice dimension. Empty
coordinate inventories correlate with empty matrix tuples.

The result does not densify matrices or recompute the Fourier sum.

## Numerical validation behavior

None. This result record retains sparse evaluator output but does not independently
recompute or numerically validate the Bloch sum.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
