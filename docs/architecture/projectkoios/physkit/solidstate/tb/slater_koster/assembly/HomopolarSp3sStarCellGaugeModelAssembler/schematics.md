# `HomopolarSp3sStarCellGaugeModelAssembler` schematics

```mermaid
flowchart TD
    I1[Validated model definition]
    I2[Two-center constructor]
    VALIDATE[Exact HomopolarSp3sStarCellGaugeModelAssembler validation]
    VALUE[Dense cell-gauge assembly result]
    I1 --> VALIDATE
    I2 --> VALIDATE
    VALIDATE --> VALUE
```

Embeds directed five-orbital hopping representatives and generated Hermitian reverses into sorted dense ordinary-cell coefficients. The action resolves groups exactly, accumulates into complex128 matrices, raises on nonfinite accumulation, and sorts displacement terms lexicographically.

The diagram shows software ownership and validation flow only. It does not imply a
research workflow, material-selection process, or evidence-promotion path.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
