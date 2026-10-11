# Module `projectkoios.physkit.solidstate.visualization.bloch_states`

## Current responsibility

This module owns unit-aware Matplotlib rendering and discrete reciprocal-index
adaptation for a sampled one-dimensional Bloch state. It distinguishes complex
amplitude, probability density, and unwrapped phase, and removes the Bloch
carrier before calculating normalized cell-periodic Fourier weights.

It does not solve an eigenproblem, normalize a state, select an eigenvector from
a degenerate subspace, infer band connectivity, or add optional widget runtime
dependencies.

## Public contracts

- [`BlochState1DPlotData`](BlochState1DPlotData/index.md)
- [`BlochState1DVisualizer`](BlochState1DVisualizer/index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/visualization/bloch_states.py`
- Tests: `tests/projectkoios/physkit/solidstate/visualization/bloch_states/`
- Laboratory: `notebooks/solidstate/bloch1d/freeelectron/free-electron-bloch-bands.ipynb`

## Navigation

- [Visualization package](../index.md)
- [Numeric contract](../numeric.md)
- [Scientific contract](../scientific.md)
