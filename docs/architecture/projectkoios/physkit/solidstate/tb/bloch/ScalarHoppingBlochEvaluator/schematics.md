# `ScalarHoppingBlochEvaluator` schematics

```mermaid
flowchart LR
    M[ScalarHoppingModel]
    Q[Unitless fractional wave vectors]
    A[ScalarHoppingBlochEvaluator]
    N[FiniteFourierSeriesEvaluator]
    V[ComplexVectorQuantity]
    R[ScalarHoppingBlochEvaluation]

    M --> A
    Q --> A
    A --> N --> V --> R
    M --> R
    Q --> R
```

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
