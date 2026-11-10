# `HomopolarSp3sStarTwoCenterIntegrals` schematics

```mermaid
flowchart TD
    I1[Seven energy quantities]
    I2[Physical unit system]
    VALIDATE[Exact HomopolarSp3sStarTwoCenterIntegrals validation]
    VALUE[Validated integral vector]
    I1 --> VALIDATE
    I2 --> VALIDATE
    VALIDATE --> VALUE
```

Retains seven real physical energy quantities in stable homopolar sp3s-star channel order. Every scalar must carry a physical energy unit, all units must match exactly, and that unit must equal the selected system energy unit.

The diagram shows software ownership and validation flow only. It does not imply a
research workflow, material-selection process, or evidence-promotion path.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
