# Package `projectkoios.physkit.periodic.lattice.embedding`

## Current responsibility

This package owns an explicit adapter from a full-column-rank ambient `3×2`
direct basis to the existing intrinsic `2×2` lattice-metric contract. It
preserves the original ambient basis, its unit, an ordered orthonormal in-plane
frame, the tangent projector, unit-aware induced metrics, and the ambient dual
reciprocal basis.

The represented object is a flat two-dimensional periodic domain with arbitrary
orientation in three-dimensional Cartesian space. It is not a full `3×3` slab
or supercell lattice.

## Public contracts

- [`EmbeddedDirectLattice2D`](EmbeddedDirectLattice2D/index.md) owns the
  immutable ambient `3×2` basis and validates length units and scale-invariant
  full column rank.
- [`EmbeddedDirectLattice2DAdaptationRequest`](EmbeddedDirectLattice2DAdaptationRequest/index.md)
  binds one nominal embedded lattice to the adapter action.
- [`EmbeddedDirectLattice2DAdapter`](EmbeddedDirectLattice2DAdapter/index.md)
  constructs the explicit ordered frame and intrinsic representation.
- [`EmbeddedDirectLattice2DAdaptation`](EmbeddedDirectLattice2DAdaptation/index.md)
  correlates ambient and intrinsic geometry and supplies explicit coordinate
  maps.
- [`EmbeddedDirectLattice2DJsonCodec`](serialization/EmbeddedDirectLattice2DJsonCodec/index.md)
  in the adjacent `embedding.serialization` module owns the strict
  schema-version-1 JSON boundary for the defining ambient basis and unit.

## Architecture views

- [Schematic](schematic/index.md)
- [Numerical](numerical/index.md)
- [Scientific](scientific/index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/embedding/__init__.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/embedding/`

## Evidence boundary

Mapped tests establish the declared representation identities, coordinate maps,
rotation/reflection invariance, finite-difference equivalence, units, and
failure behavior. They do not validate a particular material slab, atomic
geometry, or electronic model.
