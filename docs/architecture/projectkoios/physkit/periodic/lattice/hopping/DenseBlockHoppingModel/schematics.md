# `DenseBlockHoppingModel` schematics

```mermaid
flowchart TD
    ID[Model, basis, and energy identities]
    DIM[LatticeDimension]
    TERMS[Ordered DenseBlockHoppingTerm tuple]
    BASE[BaseBlockHoppingModel validation]
    MODEL[DenseBlockHoppingModel]
    ID --> MODEL
    DIM --> MODEL
    TERMS --> MODEL
    MODEL --> BASE
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
