# Package `projectkoios.physkit.periodic.lattice.rhombic`

## Current responsibility

The package owns reusable two-dimensional rhombic deformation geometry.

- [`path`](path/index.md) owns the unit-area triangular-to-square path and correlated
  point representation.
- The package façade exports reviewed path and point identities.
- No electronic model, material scale, or campaign state is attached.

## Ownership boundary

Crystallographic geometry remains separate from scientific use.

- PhysKit owns the path equation, endpoint convention, direct-basis orientation,
  metric, reciprocal geometry, and area invariant.
- Consumers own sampling, response models, thresholds, diagnostics, artifacts, and
  acceptance.
- No dependency on KSD or a calculator backend is permitted.

## Local mapping

Source, tests, and provenance mirror the rhombic package.

- Code: `src/python/projectkoios/physkit/periodic/lattice/rhombic/`
- Tests: `tests/projectkoios/physkit/periodic/lattice/rhombic/`
- [Parent lattice package](../index.md)
- [Extraction provenance](../../../../../../provenance/ksdft2effmass-crystallographic-paths.md)
