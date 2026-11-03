# Block-hopping Bloch schematics

## Representation-specific evaluation

```mermaid
flowchart LR
    DENSE[DenseBlockHoppingModel]
    SPARSE[SparseBlockHoppingModel]
    Q[Unitless fractional wave vectors]
    DE[DenseBlockHoppingBlochEvaluator]
    SE[SparseBlockHoppingBlochEvaluator]
    DR[DenseBlockHoppingBlochEvaluation]
    SR[SparseBlockHoppingBlochEvaluation]
    DM[ComplexMatrixQuantity tuple]
    SM[ComplexSparseMatrixQuantity tuple]

    DENSE --> DE
    SPARSE --> SE
    Q --> DE
    Q --> SE
    DE --> DR --> DM
    SE --> SR --> SM
```

Neither route silently crosses the dense–sparse boundary.

## Navigation

- [Module index](index.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)
