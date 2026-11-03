# `BlockHoppingTerm` numerics

`create(...)` dispatches only on exact nominal coefficient ownership. A
`ComplexMatrixQuantity` produces a dense term and a
`ComplexSparseMatrixQuantity` produces a sparse term. Unsupported quantities are
rejected.

The factory does not copy between representations, inspect numerical density, apply
a sparsity threshold, or alter values.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
