# `BaseBlockHoppingTerm` schematics

```mermaid
classDiagram
    class BaseBlockHoppingTerm {
        <<abstract>>
        +displacement LatticeDisplacement
        +coefficient_shape tuple
        +basis_size int
        +energy_unit ModelSystemUnit
        +storage_kind str
        #_validate_block_term()
    }
    class DenseBlockHoppingTerm
    class SparseBlockHoppingTerm
    BaseBlockHoppingTerm <|-- DenseBlockHoppingTerm
    BaseBlockHoppingTerm <|-- SparseBlockHoppingTerm
```

The base owns common behavior; concrete classes own exact coefficient
representations.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
