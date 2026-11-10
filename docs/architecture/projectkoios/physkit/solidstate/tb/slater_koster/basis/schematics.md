# Slater–Koster basis schematics

```mermaid
flowchart TD
    CID[Center identifier]
    TAU[Fractional center tau_c]
    CENTER[SlaterKosterBasisCenter]
    LAMBDA[Real-cubic component lambda]
    RADIAL[Radial channel r]
    COORD[Spatial coordinate c, lambda, r]
    DEF[Ordered spatial-basis definition]
    ADMIT{Exact sp3s* block order?}
    BASIS[HomopolarSp3sStarSpatialBasis]
    ROWS[Orbital-center matrix]

    CID --> CENTER
    TAU --> CENTER
    CID --> COORD
    LAMBDA --> COORD
    RADIAL --> COORD
    CENTER --> DEF
    COORD --> DEF
    DEF --> ADMIT
    ADMIT -->|yes| BASIS
    ADMIT -->|no| REJECT[Reject]
    DEF --> ROWS
```

Center order and coordinate order are observable parts of the represented state space.
The admission step compares the complete ordered coordinate tuple; it does not infer an
$sp^3s^*$ basis from orbital count, display text, or matrix dimension. The orbital-center
matrix repeats each center coordinate for every ordered orbital attached to that center.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
