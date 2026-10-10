# Embedded 2D lattice schematics

## Ownership and adaptation

```mermaid
flowchart TD
    BASIS[MatrixQuantity<br/>ambient 3 by 2 basis A]
    EMBEDDED[EmbeddedDirectLattice2D]
    REQUEST[EmbeddedDirectLattice2DAdaptationRequest]
    ADAPTER[EmbeddedDirectLattice2DAdapter]
    FRAME[Orthonormal frame Q<br/>3 by 2]
    PROJECTOR[Ambient projector P equals Q Q transpose]
    METRICS[Induced G and inverse G<br/>2 by 2]
    DUAL[Ambient dual B equals 2 pi A inverse G<br/>3 by 2]
    INTRINSIC[Intrinsic basis C equals Q transpose A<br/>2 by 2]
    METRIC[DirectLatticeMetric<br/>unchanged intrinsic contract]
    RESULT[EmbeddedDirectLattice2DAdaptation]

    BASIS --> EMBEDDED
    EMBEDDED --> REQUEST
    REQUEST --> ADAPTER
    ADAPTER --> FRAME
    FRAME --> PROJECTOR
    ADAPTER --> METRICS
    ADAPTER --> DUAL
    FRAME --> INTRINSIC
    INTRINSIC --> METRIC
    REQUEST --> RESULT
    FRAME --> RESULT
    PROJECTOR --> RESULT
    METRICS --> RESULT
    DUAL --> RESULT
    INTRINSIC --> RESULT
    METRIC --> RESULT
```

## Coordinate maps

```mermaid
flowchart LR
    FRACTIONAL[Fractional s<br/>2 components]
    AMBIENT[Ambient r equals A s<br/>3 components]
    PROJECTED_FRACTIONAL[Projected fractional<br/>r transpose B divided by 2 pi]
    INTRINSIC[Intrinsic Cartesian u<br/>2 components]
    AMBIENT_FRAME[Ambient tangent vector Q u<br/>3 components]
    PROJECTED_INTRINSIC[Projected intrinsic<br/>r transpose Q]

    FRACTIONAL -->|multiply by A| AMBIENT
    AMBIENT -->|multiply by B divided by 2 pi| PROJECTED_FRACTIONAL
    INTRINSIC -->|multiply by Q| AMBIENT_FRAME
    AMBIENT_FRAME -->|multiply by Q transpose| PROJECTED_INTRINSIC
```

## Reuse by the finite-difference capability

```mermaid
sequenceDiagram
    autonumber
    participant C as Caller
    participant E as EmbeddedDirectLattice2DAdapter
    participant R as EmbeddedDirectLattice2DAdaptation
    participant G as PeriodicFiniteDifferenceGrid
    participant L as BlochPeriodicLaplacianConstructor

    C->>E: ambient 3 by 2 basis and unit
    E->>E: construct Q, C, G, inverse G, B, P
    E-->>R: preserve ambient and intrinsic geometry
    C->>R: select intrinsic_metric and basis_unit
    C->>G: intrinsic_metric, 2D domain, basis_unit
    G-->>C: existing intrinsic 2D grid contract
    C->>L: grid and 2D Bloch twist
    L-->>C: existing quotient-seam geometric Laplacian
```

## Slab distinction

```mermaid
flowchart TD
    EMBED[Embedded 2D periodic plane<br/>A is 3 by 2]
    TWO[Two fractional periodic coordinates]
    PLANE[Rank two tangent projector]
    SLAB[Full 3D slab or supercell<br/>A is 3 by 3]
    THREE[Three fractional periodic coordinates]
    FULL[Rank three Cartesian cell]

    EMBED --> TWO
    EMBED --> PLANE
    SLAB --> THREE
    SLAB --> FULL
    EMBED -. not silently promoted to .-> SLAB
```

The adapter retains ambient orientation through `A`, `Q`, `P`, and `B`; only
the induced geometric operator is routed through the existing intrinsic `2×2`
metric contract.

## Strict serialization boundary

```mermaid
flowchart LR
    RECORD[EmbeddedDirectLattice2D]
    ENCODE[EmbeddedDirectLattice2DJsonCodec.dumps]
    WIRE[Schema v1 JSON<br/>A columns and unit only]
    DECODE[EmbeddedDirectLattice2DJsonCodec.loads]
    RESTORED[Reconstructed exact record]
    ADAPT[EmbeddedDirectLattice2DAdapter]
    DERIVED[Recomputed G, inverse G, B, Q, P]

    RECORD --> ENCODE
    ENCODE --> WIRE
    WIRE --> DECODE
    DECODE --> RESTORED
    RESTORED --> ADAPT
    ADAPT --> DERIVED
    WIRE -. excludes derived arrays .-> DERIVED
```

## Navigation

- [Module index](../index.md)
- [Numerical page](../numerical/index.md)
- [Scientific page](../scientific/index.md)
