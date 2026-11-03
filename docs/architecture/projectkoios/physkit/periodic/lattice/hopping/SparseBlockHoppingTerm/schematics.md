# `SparseBlockHoppingTerm` schematics

```mermaid
flowchart LR
    R[LatticeDisplacement R]
    Q[ComplexSparseMatrixQuantity]
    T[SparseBlockHoppingTerm]
    B[BaseBlockHoppingTerm invariants]
    R --> T
    Q --> T
    T --> B
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
