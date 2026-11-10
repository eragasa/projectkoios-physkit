# Slater–Koster bond schematics

```mermaid
flowchart TD
    R[Integer ordinary-cell displacement R]
    S[Source center identifier]
    T[Target center identifier]
    DB[SlaterKosterDirectedBond]
    DR[Physical Cartesian displacement Delta r]
    US[Physical UnitSystem]
    GEO[HomopolarSp3sStarBondGeometry]
    NORM[Euclidean normalization]
    DIR[SlaterKosterBondDirection d]
    PG[Parameter-group identifier]
    PB[Parameterized bond]
    REV[Reversed directed bond]

    R --> DB
    S --> DB
    T --> DB
    DB --> GEO
    DR --> GEO
    US --> GEO
    DR --> NORM --> DIR
    GEO --> PB
    PG --> PB
    DB --> REV
```

The ordinary displacement $R$ and physical vector $\Delta r$ remain separate inputs.
Only $R$ labels the ordinary-cell hopping coefficient. Only $\Delta r$ determines the
normalized angular direction. Reversal is a directed-bond operation; it does not mutate
or reinterpret the original record.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
