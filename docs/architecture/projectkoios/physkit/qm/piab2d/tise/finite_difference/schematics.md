# Oblique PIAB2D finite-difference schematic

```mermaid
flowchart LR
    B[MatrixQuantity direct basis] --> M[ObliquePiab2D]
    U[Mass and UnitSystem] --> M
    M --> R[ObliquePiab2DTiseFdSolveRequest]
    G[Full-grid point counts] --> R
    K[Requested low-state count] --> R
    R --> S[ObliquePiab2DTiseFdSolver]
    S --> L[Inverse-metric interior Laplacian]
    L --> H[Sparse real Hamiltonian]
    H --> E[Partial sparse or explicit dense eigensolve]
    E --> O[ObliquePiab2DTiseFdSolution]
    S --> C[Fractional and Cartesian interior coordinates]
    S --> W[Physical area weight]
    C --> O
    W --> O
```

## Ownership boundary

`ObliquePiab2D` owns bounded direct-basis geometry and coherent physical
parameters. The solve request owns one discretization choice. The solver owns
representation and diagonalization. The solution correlates every numerical
array with that model and request.

The periodic finite-difference package is not reused because its
endpoint-excluded quotient topology and seam couplings represent a different
state space. Both packages reuse direct-lattice metric geometry only.
