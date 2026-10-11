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

## Navigation

- [Package index](index.md)
- [Numeric contract](numeric.md)
- [Scientific contract](scientific.md)
