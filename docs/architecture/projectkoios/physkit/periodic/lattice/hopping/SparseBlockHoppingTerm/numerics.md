# `SparseBlockHoppingTerm` numerics

The coefficient must be an exact immutable `ComplexSparseMatrixQuantity` with
canonical CSR storage. Its shape is reported through `coefficient_shape`; the
abstract owner validates that it is nonempty and square and derives `basis_size`.

No dense materialization, tolerance-based pruning, or implicit representation
conversion occurs.

## Numerical validation behavior

None. This record retains a sparse coefficient but does not evaluate or compare it.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
