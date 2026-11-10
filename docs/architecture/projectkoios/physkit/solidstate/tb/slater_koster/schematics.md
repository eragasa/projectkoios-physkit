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
    E[Ordered basis embedding]
    H[Dense ordinary-cell coefficients]
    M[Dense block-hopping model]

    C --> B
    O --> B
    R --> G
    X --> G
    X --> D
    D --> A
    V --> T
    A --> T
    B --> E
    G --> E
    T --> E
    E --> H
    H --> M
```

The integer displacement identifies the ordinary-cell translation. The Cartesian
bond vector alone determines the normalized direction; center offsets are never folded
into the integer displacement. Assembly adds the source-to-target block at the declared
displacement and its transpose at the reversed displacement. It does not add onsite
terms or evaluate Bloch phases.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
