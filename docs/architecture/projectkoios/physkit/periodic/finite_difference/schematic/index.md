# Periodic finite-difference schematics

## Component relationships

```mermaid
flowchart TD
    DL[DirectLattice2D or DirectLattice3D]
    METRIC[DirectLatticeMetric]
    DOMAIN[FinitePeriodicDomain<br/>2D or 3D extents and C ordering]
    UNIT[PhysicalUnit or Unitless]
    GRID[PeriodicFiniteDifferenceGrid]
    LIFT[BoundaryTwistLift]
    REQUEST[BlochPeriodicLaplacianConstructionRequest]
    ACTION[BlochPeriodicLaplacianConstructor]
    REDUCTION[BoundaryTwistReductionResult]
    FIBER[TwistFiber<br/>QUOTIENT_SEAM]
    COO[Sparse COO contributions]
    CSR[ComplexSparseMatrixQuantity<br/>canonical immutable CSR]
    RESULT[BlochPeriodicLaplacian]

    DL --> METRIC
    METRIC --> GRID
    DOMAIN --> GRID
    UNIT --> GRID
    GRID --> REQUEST
    LIFT --> REQUEST
    REQUEST --> ACTION
    ACTION --> REDUCTION
    REDUCTION --> FIBER
    ACTION --> COO
    COO --> CSR
    REQUEST --> RESULT
    FIBER --> RESULT
    CSR --> RESULT
```

## Construction sequence

```mermaid
sequenceDiagram
    autonumber
    participant C as Caller
    participant G as PeriodicFiniteDifferenceGrid
    participant R as ConstructionRequest
    participant A as LaplacianConstructor
    participant T as BoundaryTwistReducer
    participant S as SciPy sparse assembly
    participant Q as ComplexSparseMatrixQuantity

    C->>G: metric, finite domain, basis unit
    G-->>C: validated endpoint-excluded grid
    C->>R: grid and BoundaryTwistLift
    R-->>C: dimension-compatible request
    C->>A: action(request=request)
    A->>T: reduce unreduced twist lift
    T-->>A: representative in [0,1) and quotient
    loop every C-order source and stencil displacement
        A->>A: divmod raw target into wrapped site and crossing q
        A->>A: coefficient times exp(2 pi i kappa dot q)
        A->>S: append COO row, column, value
    end
    A->>Q: canonicalize COO as complex CSR
    Q-->>A: immutable sparse matrix quantity
    A-->>C: BlochPeriodicLaplacian result
```

## Index and seam mechanics

```mermaid
flowchart LR
    RAW[Raw target<br/>n + delta]
    DIVMOD[Axiswise Euclidean divmod]
    WRAPPED[Wrapped coordinate<br/>inside 0 to N-1]
    QUOTIENT[Boundary quotient q]
    INDEX[C-order linear column]
    PHASE[Phase product<br/>product exp 2 pi i kappa_a q_a]
    ENTRY[COO matrix contribution]

    RAW --> DIVMOD
    DIVMOD --> WRAPPED
    DIVMOD --> QUOTIENT
    WRAPPED --> INDEX
    QUOTIENT --> PHASE
    INDEX --> ENTRY
    PHASE --> ENTRY
```

For shape `(N_0,...,N_{d-1})`, the C-order linear index is

$$
I(n_0,\ldots,n_{d-1})
=\sum_{a=0}^{d-1}n_a\prod_{b=a+1}^{d-1}N_b.
$$

A raw target is decomposed componentwise as

$$
\mathbf n+\boldsymbol\delta=\bar{\mathbf n}+N\mathbf q.
$$

The quotient-seam factor is

$$
\exp(2\pi i\,\boldsymbol\kappa\cdot\mathbf q)
=\prod_a\exp(2\pi i\kappa_aq_a).
$$

Simultaneous two-axis wraps therefore receive the product of both axis phases.

## Ownership boundary

```mermaid
flowchart LR
    GEOMETRY[periodic.finite_difference<br/>geometric Laplacian]
    ISO[Optional isotropic kinetic scaling<br/>owned by consuming quantum model]
    MASS[Optional effective-mass constitutive operator<br/>separate tensor-aware capability]
    POTENTIAL[Optional potential and Hamiltonian composition<br/>owned by consuming model]

    GEOMETRY --> ISO
    GEOMETRY -. not a scalar matrix multiplication .-> MASS
    ISO --> POTENTIAL
```

The package ends at the geometric Laplacian. Kinetic scaling, anisotropic
constitutive tensors, potentials, and material interpretation remain separate.

## Navigation

- [Package index](../index.md)
- [Numerical page](../numerical/index.md)
- [Scientific page](../scientific/index.md)
