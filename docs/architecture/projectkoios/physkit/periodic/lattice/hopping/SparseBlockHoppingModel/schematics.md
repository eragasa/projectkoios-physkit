# `SparseBlockHoppingModel` schematics

```mermaid
flowchart TD
    ID[Model, basis, and energy identities]
    DIM[LatticeDimension]
    TERMS[Ordered SparseBlockHoppingTerm tuple]
    BASE[BaseBlockHoppingModel validation]
    MODEL[SparseBlockHoppingModel]
    ID --> MODEL
    DIM --> MODEL
    TERMS --> MODEL
    MODEL --> BASE
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
