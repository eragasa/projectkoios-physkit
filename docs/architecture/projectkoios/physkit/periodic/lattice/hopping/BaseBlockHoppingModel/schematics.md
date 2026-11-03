# `BaseBlockHoppingModel` schematics

```mermaid
classDiagram
    class BaseBlockHoppingModel {
        <<abstract>>
        +identifier str
        +dimension LatticeDimension
        +terms tuple
        +energy_unit
        +energy_reference str
        +basis_identifier str
        +basis_size int
        +storage_kind str
        #_validate_block_model()
    }
    class DenseBlockHoppingModel
    class SparseBlockHoppingModel
    BaseBlockHoppingModel <|-- DenseBlockHoppingModel
    BaseBlockHoppingModel <|-- SparseBlockHoppingModel
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
