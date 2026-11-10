# Slater–Koster parameter-group schematics

```mermaid
flowchart LR
    ID[Nonempty group identifier]
    V[Seven validated energy integrals]
    GROUP[HomopolarSp3sStarParameterGroup]
    B1[Parameterized bond 1]
    B2[Parameterized bond 2]
    BN[Parameterized bond n]

    ID --> GROUP
    V --> GROUP
    GROUP -. explicit sharing .-> B1
    GROUP -. explicit sharing .-> B2
    GROUP -. explicit sharing .-> BN
```

The dotted edges denote declared reuse of one immutable integral vector. They do not
mean that the bonds share a geometric shell, length, chemical environment, or material
interpretation unless a client explicitly establishes those facts.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
