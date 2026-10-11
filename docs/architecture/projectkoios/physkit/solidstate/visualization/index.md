# Package `projectkoios.physkit.solidstate.visualization`

This package owns reusable Matplotlib rendering behavior for solid-state
results and plot-data records. Visualizers accept or create `Axes`, draw without
calling `show`, and return the owning axes so notebooks and document scripts
retain layout control.

`BandStructurePlotData` represents any sampled band energies along an ordered
scalar path coordinate. The coordinate may parameterize a one-, two-, or
three-dimensional reciprocal-space path; the visualizer does not interpret or
reconstruct reciprocal vectors. `BandStructureVisualizer` plots every energy
column and can add a same-grid reference matrix after checking shape and energy
units.

Model-specific code owns conversion from a physical result to the scalar path
coordinate. High-symmetry labels, path construction, material interpretation,
and figure-specific annotations remain with the caller unless a reusable typed
contract is introduced.

Older electronic- and phonon-specific visualizers remain supported by their
defining modules. The dimension-independent band visualizer does not change
their contracts.

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/visualization/`
- Tests: `tests/projectkoios/physkit/solidstate/visualization/`

## Navigation

- [Schematics](schematics.md)
- [Numeric contract](numeric.md)
- [Scientific contract](scientific.md)
