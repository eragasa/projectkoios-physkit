# `ScalarHoppingBlochEvaluation` schematics

```mermaid
flowchart LR
    M[Exact ScalarHoppingModel]
    Q[Exact fractional-wave-vector MatrixQuantity]
    V[ComplexVectorQuantity]
    R[ScalarHoppingBlochEvaluation]
    M --> R
    Q --> R
    V --> R
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
