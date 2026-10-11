# Density-of-states schematics

## Ownership and evaluation

```mermaid
flowchart LR
    MESH[Unique equal-weight<br/>full-zone samples]
    BANDS[MatrixQuantity<br/>sample x band energies]
    MODEL[dos1d, dos2d, or dos3d<br/>BZ measure and degeneracy]
    GRID[Energy grid or bin edges]
    GAUSS[Gaussian evaluator]
    HIST[Histogram evaluator]
    RESULT[Unit-aware DOS result<br/>represented and integrated states]
    VIEW[Dimension-independent<br/>DOS visualizer]

    MESH --> BANDS
    BANDS --> GAUSS
    BANDS --> HIST
    MODEL --> GAUSS
    MODEL --> HIST
    GRID --> GAUSS
    GRID --> HIST
    GAUSS --> RESULT
    HIST --> RESULT
    RESULT --> VIEW
```

Mesh generation and band solving remain with the dimensional physical model.
The DOS owner starts from the resulting sampled energies.

## Dimensional normalization

```mermaid
flowchart TD
    D1[dos1d<br/>BZ reciprocal length]
    D2[dos2d<br/>BZ reciprocal area]
    D3[dos3d<br/>BZ reciprocal volume]
    SHARED[Shared equal-weight state factor<br/>g / ((2 pi)^d N)]
    OUT[DOS per energy per<br/>length, area, or volume]

    D1 --> SHARED
    D2 --> SHARED
    D3 --> SHARED
    SHARED --> OUT
```

## Navigation

- [Package index](index.md)
- [Numeric contract](numeric.md)
- [Scientific contract](scientific.md)
