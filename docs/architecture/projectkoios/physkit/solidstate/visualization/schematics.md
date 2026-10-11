# Solid-state visualization schematics

## Dimension-independent band rendering

```mermaid
flowchart LR
    R1[1D band result]
    R2[2D reciprocal-space path result]
    R3[3D reciprocal-space path result]
    ADAPT[Model-owned path adaptation]
    PATH[VectorQuantity scalar path coordinate]
    ENERGY[MatrixQuantity band energies]
    DATA[BandStructurePlotData]
    VIEW[BandStructureVisualizer]
    AX[Matplotlib Axes]

    R1 --> ADAPT
    R2 --> ADAPT
    R3 --> ADAPT
    ADAPT --> PATH
    ADAPT --> ENERGY
    PATH --> DATA
    ENERGY --> DATA
    DATA --> VIEW
    VIEW --> AX
```

The visualizer consumes a scalar parameterization of an already selected path.
It does not inspect the dimension of the reciprocal vectors used to construct
that path.

## Optional reference overlay

```mermaid
flowchart LR
    PRIMARY[Primary band-energy matrix]
    REFERENCE[Optional reference matrix]
    CHECK[Shape and energy-unit checks]
    DRAW[Primary solid lines<br/>reference dashed lines]

    PRIMARY --> CHECK
    REFERENCE --> CHECK
    CHECK --> DRAW
```

Notebook or model code owns the scientific meaning of the reference, such as an
exact free-particle spectrum or a convergence baseline.

## One-dimensional sampled Bloch states

```mermaid
flowchart LR
    POS[Uniform endpoint-excluded positions]
    PSI[Complex Bloch amplitudes psi_k]
    K[Explicit Bloch wave number k]
    DATA[BlochState1DPlotData]
    STATE[Amplitude, density, phase]
    REMOVE[Remove exp(ikx)]
    FFT[Cell-periodic FFT weights]

    POS --> DATA
    PSI --> DATA
    K --> DATA
    DATA --> STATE
    DATA --> REMOVE --> FFT
```

The full twisted state is rendered in real space. Reciprocal integer labels
belong to the periodic factor $u_k=e^{-ikx}\psi_k$, not directly to the
nonperiodic carrier-bearing samples.

## Navigation

- [Package index](index.md)
- [Numeric contract](numeric.md)
- [Scientific contract](scientific.md)
