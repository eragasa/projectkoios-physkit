# Slater–Koster assembly schematics

```mermaid
flowchart TD
    BASIS[Exact ordered sp3s* basis]
    GROUPS[Unique parameter groups]
    BONDS[Unique directed bond representatives]
    DEF[Cell-gauge model definition]
    LOOKUP[Resolve group by exact identifier]
    T[Construct 5 x 5 hopping T]
    EMBED[Embed source-to-target block]
    REVERSE[Generate transpose reverse block]
    ACCUM[Accumulate by integer displacement]
    SORT[Lexicographic displacement sort]
    TERMS[DenseBlockHoppingTerm tuple]
    MODEL[DenseBlockHoppingModel]
    RESULT[Correlated assembly result]

    BASIS --> DEF
    GROUPS --> DEF
    BONDS --> DEF
    DEF --> LOOKUP --> T
    BASIS --> EMBED
    BONDS --> EMBED
    T --> EMBED
    EMBED --> ACCUM
    T --> REVERSE --> ACCUM
    ACCUM --> SORT --> TERMS --> MODEL --> RESULT
    DEF --> RESULT
```

Forward and reverse contributions may target the same accumulator when the ordinary
cell displacement is zero and the centers differ. Accumulation occurs before immutable
quantity construction. The generated model retains basis and energy-reference metadata,
but the structured definition remains the authoritative exact basis correlation.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
