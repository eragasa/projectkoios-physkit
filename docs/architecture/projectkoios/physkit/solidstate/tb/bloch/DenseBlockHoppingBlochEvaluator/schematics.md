# `DenseBlockHoppingBlochEvaluator` schematics

```mermaid
flowchart TD
    M[DenseBlockHoppingModel]
    Q[MatrixQuantity fractional wave vectors]
    MODE{Configured backend}
    ORDERED[Entrywise ordered finite series]
    BATCHED[Shared-grid batched finite series]
    VALUES[Dense complex matrix stack]
    RESULT[DenseBlockHoppingBlochEvaluation]

    M --> MODE
    Q --> MODE
    MODE --> ORDERED --> VALUES
    MODE --> BATCHED --> VALUES
    VALUES --> RESULT
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
