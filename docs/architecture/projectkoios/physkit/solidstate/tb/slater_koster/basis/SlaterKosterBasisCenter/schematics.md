# `SlaterKosterBasisCenter` schematics

```mermaid
flowchart TD
    I1[Center identifier]
    I2[Fractional coordinate]
    VALIDATE[Exact SlaterKosterBasisCenter validation]
    VALUE[Basis center]
    I1 --> VALIDATE
    I2 --> VALIDATE
    VALIDATE --> VALUE
```

Binds one nonempty center identifier to an immutable primitive fractional coordinate. The coordinate shape is exactly `(3,)`; no Cartesian lattice basis or unit conversion is accepted at this boundary.

The diagram shows software ownership and validation flow only. It does not imply a
research workflow, material-selection process, or evidence-promotion path.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
