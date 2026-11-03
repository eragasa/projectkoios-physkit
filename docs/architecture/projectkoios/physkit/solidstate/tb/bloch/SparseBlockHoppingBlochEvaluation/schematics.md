# `SparseBlockHoppingBlochEvaluation` schematics

```mermaid
flowchart LR
    M[Exact SparseBlockHoppingModel]
    Q[Exact fractional-wave-vector MatrixQuantity]
    H[Ordered ComplexSparseMatrixQuantity tuple]
    R[SparseBlockHoppingBlochEvaluation]
    M --> R
    Q --> R
    H --> R
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
