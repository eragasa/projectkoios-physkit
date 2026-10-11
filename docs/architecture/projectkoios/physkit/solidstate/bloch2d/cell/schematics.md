# `BlochCell2D` schematics

## Ownership and derived geometry

```mermaid
flowchart TD
    BASIS[MatrixQuantity direct_basis<br/>ordered 2 x 2 columns]
    SYSTEM[UnitSystem.METAL or<br/>UnitSystem.NONDIMENSIONAL]
    CELL[BlochCell2D]
    RANK[Scale-invariant<br/>normalized-area check]
    AREA[ScalarQuantity area]
    RECIP[MatrixQuantity reciprocal_basis]
    METRIC[MatrixQuantity metric_tensor]
    CHECKER[CheckerboardPotentialBloch2D]

    BASIS --> CELL
    SYSTEM --> CELL
    CELL --> RANK
    RANK --> AREA
    RANK --> RECIP
    RANK --> METRIC
    CELL --> CHECKER
```

`BlochCell2D` owns represented geometry and units. It does not own a potential,
particle mass convention, reciprocal cutoff, band filling, or material label.

## Construction and validation sequence

```mermaid
sequenceDiagram
    autonumber
    participant C as Caller
    participant B as BlochCell2D
    participant V as normalized-column-volume utility
    participant Q as Derived quantities

    C->>B: direct_basis and unit_system
    B->>B: require MatrixQuantity and UnitSystem
    B->>B: allow METAL or NONDIMENSIONAL only
    B->>B: match basis unit to selected system
    B->>B: require shape (2, 2)
    B->>V: test normalized column area > 1e-8
    V-->>B: rank decision independent of uniform scale
    B-->>C: immutable validated cell
    C->>Q: request area, reciprocal basis, or metric
    Q-->>C: quantity in squared, reciprocal, or squared basis units
```

## Consumer boundary

```mermaid
flowchart LR
    CELL[BlochCell2D]
    SHAPE[Scale-normalized square-cell check]
    MODEL[CheckerboardPotentialBloch2D]
    METAL[METAL convention<br/>angstrom, eV, dalton]
    ND[NONDIMENSIONAL convention<br/>hbar = me = 1]
    CSR[Canonical complex CSR Hamiltonian]
    DETECT[Sparse-operator detection]
    ARPACK[Sparse partial eigensolver]
    DENSE[Explicit dense fallback<br/>complete or all-but-one spectrum]

    CELL --> SHAPE
    SHAPE --> MODEL
    METAL --> MODEL
    ND --> MODEL
    MODEL --> CSR
    CSR --> DETECT
    DETECT --> ARPACK
    DETECT --> DENSE
```

The checkerboard model adds a square-shape restriction. General `BlochCell2D`
geometry may be skew and is not itself restricted to square cells. The
Hamiltonian representation remains sparse by default; only eigensolver requests
outside the maintained ARPACK partial-spectrum range materialize a dense copy.

## Navigation

- [Module index](index.md)
- [Numeric page](numeric.md)
- [Scientific page](scientific.md)
