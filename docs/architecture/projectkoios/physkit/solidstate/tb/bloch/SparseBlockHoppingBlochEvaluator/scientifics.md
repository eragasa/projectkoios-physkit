# `SparseBlockHoppingBlochEvaluator` scientific contract

The evaluator assembles the same Bloch sum as the dense route while preserving sparse
internal storage. The output sparsity pattern is the represented phase-weighted union
of source block patterns; sparsity is not guaranteed to remain small for every model.

The action does not infer Hermiticity, diagonalize the result, or establish that a
sparse eigensolver is appropriate.

- [Index](index.md)
- [Schematics](schematics.md)
- [Numerics](numerics.md)
