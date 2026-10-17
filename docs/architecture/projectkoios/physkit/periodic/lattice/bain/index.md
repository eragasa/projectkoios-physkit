# Package `projectkoios.physkit.periodic.lattice.bain`

## Current responsibility

The package owns the reusable crystallographic Bain family.

- [`bridge`](bridge/index.md) owns exact coordinate conversion between declared FCC and
  Bain primitive bases.
- [`path`](path/index.md) owns the unit-volume Bravais Bain path and species-neutral
  diamond decoration.
- The package façade exports reviewed bridge, path, and point identities.
- Electronic response and material-transition interpretation remain outside the package.

## Ownership boundary

The package stops at represented crystallographic geometry.

- PhysKit owns basis conventions, coordinate maps, deformation equations, fractional
  decorations, and correlated lattice/metric evidence.
- Consumers own path sampling, material parameters, Hamiltonians, response models,
  thresholds, artifacts, and acceptance.
- No dependency on KSD or a calculator backend is permitted.

## Local mapping

Source, tests, and provenance mirror the Bain package.

- Code: `src/python/projectkoios/physkit/periodic/lattice/bain/`
- Tests: `tests/projectkoios/physkit/periodic/lattice/bain/`
- [Parent lattice package](../index.md)
- [Extraction provenance](../../../../../../provenance/ksdft2effmass-crystallographic-paths.md)
