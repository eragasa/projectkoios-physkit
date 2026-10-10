# Package `projectkoios.physkit.periodic.lattice`

## Current responsibility

The package owns physical direct-lattice records, direct-basis metric geometry,
and reusable finite-periodic integer-domain mechanics. The finite-periodic slice
includes boundary twists, integral symmetry operations, scalar hopping models,
represented sparse operators, compatibility-gated composition, uniform-link
construction, and explicit gauge bridges.

## Public contract

Physical `DirectLattice1D`, `DirectLattice2D`, and `DirectLattice3D` records are
not interchangeable with `FinitePeriodicDomain`. Supported finite-periodic
interfaces are intentionally limited to one, two, and three dimensions.
Defining modules remain the preferred import routes.

## Navigation

- [Module `lattice2d`](lattice2d/index.md)
- [Module `lattice3d`](lattice3d/index.md)
- [Module `metric`](metric/index.md)
- [Module `gauge_bridges`](gauge_bridges/index.md)
- [Parent package](../index.md)
- [Source provenance](../../../../../provenance/finite-periodic-lattice-source-mapping.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/__init__.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/`

## Evidence

Implementation conformance is supported by the direct-basis metric,
nearest-image, finite-domain, boundary-phase, symmetry, hopping,
represented-operator, composition, construction, and gauge-bridge tests. These
tests make no campaign-level or scientific-validation claim.
