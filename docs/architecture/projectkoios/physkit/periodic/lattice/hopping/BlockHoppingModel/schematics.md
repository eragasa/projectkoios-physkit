# `BlockHoppingModel` schematics

```mermaid
flowchart LR
    META[Model identities and unit]
    TERMS{Homogeneous exact term inventory}
    DENSE[DenseBlockHoppingTerm tuple]
    SPARSE[SparseBlockHoppingTerm tuple]
    DM[DenseBlockHoppingModel]
    SM[SparseBlockHoppingModel]

    META --> TERMS
    DENSE --> TERMS
    SPARSE --> TERMS
    TERMS -->|dense| DM
    TERMS -->|sparse| SM
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
