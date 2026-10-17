# Package `projectkoios.physkit.periodic.lattice`

## Current responsibility

The package owns reusable physical and represented lattice contracts.

- Physical geometry includes direct lattices, direct-basis metrics, named
  crystallographic deformation paths, and explicit ambient-3D/intrinsic-2D adaptation.
- Finite-periodic mechanics include integer domains, boundary twists, integral
  symmetry operations, scalar and block hopping models, and represented sparse
  operators.
- Composition capabilities include compatibility gates, uniform-link
  construction, and explicit gauge bridges.

## Public contract

Physical `DirectLattice1D`, `DirectLattice2D`, and `DirectLattice3D` records are
not interchangeable with `FinitePeriodicDomain`. Supported finite-periodic
interfaces are intentionally limited to one, two, and three dimensions.

- The package façade exports direct-lattice records, `DirectLatticeMetric`, and
  reviewed crystallographic path records/actions as supported composition surfaces.
- `DirectLatticeMetric` retains its defining-module identity and current 2D/3D
  semantics.
- Defining modules remain the preferred routes for the remaining interfaces.

## Navigation

The architecture pages follow the package's owned modules and subpackages.

- [Module `deformation_paths`](deformation_paths/index.md)
- [Module `lattice2d`](lattice2d/index.md)
- [Module `lattice3d`](lattice3d/index.md)
- [Package `embedding`](embedding/index.md)
- [Module `metric`](metric/index.md)
- [Module `hopping`](hopping/index.md)
- [Module `gauge_bridges`](gauge_bridges/index.md)
- [Parent package](../index.md)
- [Source provenance](../../../../../provenance/finite-periodic-lattice-source-mapping.md)

## Local mapping

Maintained implementation and verification mirror this package ownership.

- Code: `src/python/projectkoios/physkit/periodic/lattice/__init__.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/`

## Evidence

Mapped tests provide implementation-conformance evidence for the owned contracts.

- Coverage includes direct-basis metrics, crystallographic path invariants, nearest
  images, finite domains, boundary phases, symmetry, hopping, represented operators,
  composition, construction, and gauge bridges.
- These checks make no campaign-level or scientific-validation claim.
