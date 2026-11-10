# Slater–Koster channel-basis schematics

```mermaid
flowchart TD
    D[Normalized direction d = l,m,n]
    ORDER[Stable seven-channel order]
    SS[ss sigma matrix]
    SP[sp sigma matrix]
    PPS[pp sigma matrix]
    PPP[pp pi matrix]
    SSTAR[s* p sigma matrix]
    SSSTAR[s s* sigma matrix]
    SSTARSTAR[s* s* sigma matrix]
    RESULT[Ordered channel-basis result]

    D --> SS
    D --> SP
    D --> PPS
    D --> PPP
    D --> SSTAR
    D --> SSSTAR
    D --> SSTARSTAR
    ORDER --> RESULT
    SS --> RESULT
    SP --> RESULT
    PPS --> RESULT
    PPP --> RESULT
    SSTAR --> RESULT
    SSSTAR --> RESULT
    SSTARSTAR --> RESULT
```

Every branch constructs one independent coefficient matrix. The result order is the
same order used by integral records and direct hopping construction, so no string-based
channel matching or inferred shell semantics enters the linear combination.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
