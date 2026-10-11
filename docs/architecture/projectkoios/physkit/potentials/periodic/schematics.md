# Periodic-potential schematics

## Ownership boundary

```mermaid
flowchart LR
    AMP[ScalarQuantity amplitude]
    SHAPE[Harmonics and phases]
    POT[BasePeriodicPotential]
    COORD[Unitless fractional coordinates]
    MODES[Ordered integer reciprocal modes]
    VALUES[VectorQuantity potential values]
    COEFF[ComplexVectorQuantity Fourier coefficients]

    AMP --> POT
    SHAPE --> POT
    COORD --> POT
    MODES --> POT
    POT --> VALUES
    POT --> COEFF
```

The potential owns its scalar field and Fourier convention. The lattice owner
maps physical positions to fractional coordinates and supplies reciprocal-mode
identifiers.

## Potential-to-operator conversion

```mermaid
flowchart TD
    POT[Periodic potential]
    SAMPLE[Fractional sample coordinates]
    BASIS[Ordered reciprocal-mode basis]
    REAL[Sampled-diagonal converter]
    FOURIER[Fourier-operator converter]
    DIAG[Canonical real CSR diagonal operator]
    CONV[Canonical complex CSR convolution operator]
    FD[Finite-difference Hamiltonian owner]
    PW[Plane-wave Hamiltonian owner]

    POT --> REAL
    SAMPLE --> REAL
    REAL --> DIAG
    DIAG --> FD

    POT --> FOURIER
    BASIS --> FOURIER
    FOURIER --> CONV
    CONV --> PW
```

A converter does not add kinetic energy or solve an eigenproblem. Hamiltonian
owners compose the returned potential operator with their separately represented
kinetic operator.

## Current consumers

```mermaid
flowchart LR
    C1[CosinePotential1D]
    C2[CheckerboardPotential2D]
    B1[CosinePotentialBloch1D]
    B2[CheckerboardPotentialBloch2D]
    D[Sampled diagonal converter]
    F[Fourier converter]

    C1 --> D --> B1
    C2 --> F --> B2
```

The three-dimensional cosine-product model is currently a reusable field and
operator source; it does not imply that a three-dimensional Bloch solver has
been selected or validated.

## Navigation

- [Package index](index.md)
- [Numeric contract](numeric.md)
- [Scientific contract](scientific.md)
