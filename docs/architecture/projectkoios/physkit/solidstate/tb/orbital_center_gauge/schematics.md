# Orbital-center gauge schematics

```mermaid
flowchart TD
    TAU[Ordered fractional centers Tau]
    Q[Fractional wave vector q]
    DEF[Gauge transformation definition]
    DOT[Rowwise dot products Tau q]
    PHASE[Positive phases exp +2 pi i Tau q]
    DIAG[Diagonal matrix D q]
    RESULT[Correlated gauge result]
    HCELL[Cell-gauge matrix H_cell q]
    LEFT[Left multiply by D dagger]
    HFULL[Full-center matrix H_full q]

    TAU --> DEF
    Q --> DEF
    DEF --> DOT --> PHASE --> DIAG --> RESULT
    HCELL --> LEFT
    DIAG --> LEFT
    LEFT --> HFULL
    DIAG --> HFULL
```

The ordered center rows and the represented matrix basis must refer to the same exact
coordinate order. The module retains the order supplied by the caller and performs no
sorting, center lookup, reciprocal-axis permutation, or compatibility inference.

The transformation changes representation rather than the represented operator. Both
left and right factors are required; multiplying on only one side is not a basis change.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
