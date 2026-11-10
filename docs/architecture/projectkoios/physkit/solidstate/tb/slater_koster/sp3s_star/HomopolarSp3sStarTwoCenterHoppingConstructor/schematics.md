# `HomopolarSp3sStarTwoCenterHoppingConstructor` schematics

```mermaid
flowchart TD
    I1[Hopping definition]
    I2[Angular channel basis]
    VALIDATE[Exact HomopolarSp3sStarTwoCenterHoppingConstructor validation]
    VALUE[Two-center hopping result]
    I1 --> VALIDATE
    I2 --> VALIDATE
    VALIDATE --> VALUE
```

Constructs one energy-valued two-center matrix by ordered linear combination of the seven angular channel matrices. Binary64 accumulation is checked for nonfinite output. Overflow raises rather than clipping or silently regularizing values.

The diagram shows software ownership and validation flow only. It does not imply a
research workflow, material-selection process, or evidence-promotion path.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
