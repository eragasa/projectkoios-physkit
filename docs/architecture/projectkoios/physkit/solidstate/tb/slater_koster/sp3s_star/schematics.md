# Homopolar $sp^3s^*$ hopping schematics

```mermaid
flowchart TD
    D[SlaterKosterBondDirection]
    V[Seven energy integrals]
    DEF[Two-center hopping definition]
    B[Seven nondimensional channel matrices]
    SCALE[Scale each B_alpha by V_alpha]
    SUM[Binary64 matrix sum]
    FINITE{All entries finite?}
    RESULT[Two-center hopping result]
    OVERFLOW[Raise OverflowError]

    D --> DEF
    V --> DEF
    D --> B
    DEF --> SCALE
    B --> SCALE
    SCALE --> SUM --> FINITE
    FINITE -->|yes| RESULT
    FINITE -->|no| OVERFLOW
```

The constructor delegates angular placement to the channel-basis owner and performs one
ordered linear reduction. It does not infer absent channels: a missing physical
interaction is represented by an explicit zero-valued energy quantity.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
