# Package `projectkoios.physkit.solidstate.visualization`

## Current responsibility

The package owns reusable Matplotlib rendering behavior for solid-state results and
plot-data records.

- Visualizers accept or create `Axes`, draw without calling `show`, and return the
  owning axes.
- `BandStructurePlotData` represents sampled band energies along an ordered scalar path
  coordinate without reconstructing reciprocal vectors.
- The [`bloch_states`](bloch_states/index.md) family owns unit-aware rendering for
  sampled one-dimensional Bloch states.
- Older electronic- and phonon-specific visualizers remain supported by their defining
  modules.

## Composition boundary

Model-specific code prepares physical results before rendering.

- The dimension-independent band visualizer accepts one-, two-, or three-dimensional
  paths after conversion to a scalar path coordinate.
- Its optional reference matrix must share the primary grid, shape, and energy unit.
- Older 1D electronic and phonon visualizers consume the canonical
  `periodic.lattice.ReciprocalLattice1D` and normalize coordinates by its declared
  `b1` basis vector.
- High-symmetry labels, path construction, materials, and figure-specific annotations
  remain with the caller unless a reusable typed contract owns them.

## Local mapping

Maintained code and tests mirror visualization ownership.

- Code: `src/python/projectkoios/physkit/solidstate/visualization/`
- Tests: `tests/projectkoios/physkit/solidstate/visualization/`

## Navigation

Supporting pages document reusable rendering contracts.

- [Schematics](schematics.md)
- [Numeric contract](numeric.md)
- [Scientific contract](scientific.md)
