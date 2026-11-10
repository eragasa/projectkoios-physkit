# `SlaterKosterSpatialBasisCoordinate` schematics

```mermaid
flowchart TD
    I1[Center identifier]
    I2[Orbital component]
    I3[Radial channel]
    VALIDATE[Exact SlaterKosterSpatialBasisCoordinate validation]
    VALUE[Spatial basis coordinate]
    I1 --> VALIDATE
    I2 --> VALIDATE
    I3 --> VALIDATE
    VALIDATE --> VALUE
```

Identifies one ordered spin-independent spatial coordinate by center, angular component, and radial channel. Booleans, NumPy integers, negative channels, empty center labels, and arbitrary component strings are rejected.

The diagram shows software ownership and validation flow only. It does not imply a
research workflow, material-selection process, or evidence-promotion path.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
