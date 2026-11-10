# Slater–Koster schematics

```mermaid
flowchart TD
    C[Ordered centers]
    O[Orbital component and radial channel]
    B[Homopolar sp3s* spatial basis]
    R[Integer ordinary-cell displacement]
    X[Unit-bearing Cartesian bond vector]
    G[Directed bond geometry]
    D[Normalized Cartesian bond direction]
    V[Seven energy integrals]
    A[Seven nondimensional angular matrices]
    T[Real 5 x 5 two-center hopping]

    C --> B
    O --> B
    R --> G
    X --> G
    X --> D
    D --> A
    V --> T
    A --> T
```

The integer displacement identifies the ordinary-cell translation. The Cartesian
bond vector alone determines the normalized direction; center offsets are never folded into
the integer displacement.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
