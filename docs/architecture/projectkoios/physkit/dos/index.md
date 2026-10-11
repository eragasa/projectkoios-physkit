# Package `projectkoios.physkit.dos`

This package owns shared contracts and numerical estimators for densities of
states derived from equally weighted, unique samples over a complete Brillouin
zone.

`BaseUniformBandDensityOfStates` owns reciprocal-space dimension,
Brillouin-zone measure, and state degeneracy. Its Gaussian and histogram
facades consume a `(sample_count, band_count)` energy matrix and return
unit-aware density-of-states results together with represented and numerically
included state densities.

`BaseUniformBandDensityOfStates.solve_chemical_potential(...)` maps every
sampled eigenvalue to its exact full-zone state-density weight and delegates the
finite-temperature inversion to the weighted Fermi--Dirac solver. It does not
solve on histogram centers or a broadened display grid.

The [`thermal`](thermal/index.md) module transforms a normalized DOS with one
`FermiDiracState` into occupation factors, occupied and unoccupied DOS, the
thermal window, and thermally active DOS.

`DensityOfStatesPlotData` adapts ordinary or thermal results to correlated
energy and density vectors. `DensityOfStatesVisualizer` renders energy on the
horizontal or vertical axis without assuming reciprocal-space dimension and
without choosing normalization, bins, broadening, or occupations.

Concrete normalization models live in:

- [`dos1d`](../dos1d/index.md),
- [`dos2d`](../dos2d/index.md), and
- [`dos3d`](../dos3d/index.md).

The package does not construct reciprocal meshes, infer whether endpoint
samples are duplicated, choose occupied bands, assign a material, or solve a
Hamiltonian. Callers must provide an endpoint-excluded or otherwise unique
full-zone mesh when using equal sample weights.

## Local mapping

- Code: `src/python/projectkoios/physkit/dos/`
- Tests: `tests/projectkoios/physkit/dos/`

## Navigation

- [Schematics](schematics.md)
- [Numeric contract](numeric.md)
- [Scientific contract](scientific.md)
