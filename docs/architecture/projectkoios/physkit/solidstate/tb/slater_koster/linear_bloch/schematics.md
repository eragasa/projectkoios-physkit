# Linear Slater–Koster Bloch schematics

```mermaid
flowchart LR
    PARAM[Linear parameterization]
    FLAT[Flatten model family internally]
    Q[Ordered reciprocal coordinates]
    BLOCH[DenseBlockHoppingBlochEvaluator]
    REGROUP[Regroup without flattening matrix entries]
    RESULT[Group / channel / sample matrices]

    PARAM --> FLAT --> BLOCH
    Q --> BLOCH
    BLOCH --> REGROUP --> RESULT
    PARAM --> RESULT
    Q --> RESULT
```

Only the model-family axis is flattened for delegation. Matrix entries remain matrices,
and the returned result restores the exact group/channel nesting.

- [Index](index.md)
- [Numerics](numerics.md)
- [Scientific contract](scientifics.md)
