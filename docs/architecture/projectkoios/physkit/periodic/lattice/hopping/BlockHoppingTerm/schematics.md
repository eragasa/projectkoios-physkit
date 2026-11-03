# `BlockHoppingTerm` schematics

```mermaid
flowchart LR
    D[LatticeDisplacement]
    C{Exact coefficient type}
    DM[ComplexMatrixQuantity]
    SM[ComplexSparseMatrixQuantity]
    DT[DenseBlockHoppingTerm]
    ST[SparseBlockHoppingTerm]

    D --> C
    DM --> C
    SM --> C
    C -->|dense| DT
    C -->|sparse| ST
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
