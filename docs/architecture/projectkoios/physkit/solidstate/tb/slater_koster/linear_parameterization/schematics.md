# Linear Slater–Koster parameterization schematics

```mermaid
flowchart TD
    DEF[Admitted cell-gauge definition]
    GROUP[Parameter group g]
    CHANNEL[Stable channel alpha]
    BONDS[Bonds explicitly naming g]
    ANGULAR[Angular matrix B_alpha of d]
    EMBED[Embed source-target block]
    REVERSE[Generate transpose reverse]
    ACCUM[Accumulate by integer displacement]
    MODEL[Nondimensional DenseBlockHoppingModel]
    RESULT[Group-major channel-model result]

    DEF --> GROUP
    DEF --> BONDS
    GROUP --> BONDS
    CHANNEL --> ANGULAR
    BONDS --> ANGULAR --> EMBED
    EMBED --> ACCUM
    ANGULAR --> REVERSE --> ACCUM
    ACCUM --> MODEL --> RESULT
```

Parameter values are not read during construction. Group identifiers select explicit
sharing only and never infer shell or material meaning.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
