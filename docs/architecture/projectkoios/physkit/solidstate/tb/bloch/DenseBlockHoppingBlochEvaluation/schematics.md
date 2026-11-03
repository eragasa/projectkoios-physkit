# `DenseBlockHoppingBlochEvaluation` schematics

```mermaid
flowchart LR
    M[Exact DenseBlockHoppingModel]
    Q[Exact fractional-wave-vector MatrixQuantity]
    H[Ordered ComplexMatrixQuantity tuple]
    R[DenseBlockHoppingBlochEvaluation]
    M --> R
    Q --> R
    H --> R
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
