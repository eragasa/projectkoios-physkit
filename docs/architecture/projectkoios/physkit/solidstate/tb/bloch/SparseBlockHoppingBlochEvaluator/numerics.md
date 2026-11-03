# `SparseBlockHoppingBlochEvaluator` numerics

For each coordinate row, the evaluator starts with an empty complex128 CSR array and
adds every coefficient CSR block multiplied by its complex phase. The result is
canonicalized by summing duplicates, eliminating exact zeros, and sorting indices.
Nonfinite real or imaginary data raises `FloatingPointError`.

The evaluator does not construct a dense matrix, call the dense evaluator, apply a
sparsity threshold, or provide a dense fallback.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
