# Module `projectkoios.physkit.periodic.lattice.deformation_paths`

## Current responsibility

The module owns reusable crystallographic deformation-path representations and laws.

- [`UnitVolumeBainPath3D`](UnitVolumeBainPath3D/index.md) evaluates one explicit
  FCC-to-BCC Bain convention in a body-centered-tetragonal primitive basis.
- [`UnitVolumeBainPathPoint3D`](UnitVolumeBainPathPoint3D/index.md) retains correlated
  scalar and direct-lattice evidence at one path parameter.
- The package façade deliberately exports both identities.
- Later named crystallographic paths remain separate bounded contracts.

## Bain convention

The represented path fixes its mathematical convention rather than selecting a
campaign.

- The closed path parameter runs from zero at the FCC endpoint to one at the BCC
  endpoint.
- The longitudinal-to-transverse ratio varies linearly from `sqrt(2)` to one.
- Transverse and longitudinal scales preserve unit primitive-cell volume.
- Primitive vectors use one declared body-centered-tetragonal orientation and ordering.
- Direct metrics and reciprocal geometry derive from the exact represented lattice.

## Data flow

The path action returns owner-native lattice geometry without campaign state.

- Evaluation constructs one correlated immutable point.
- The point owns a validated copy of its direct lattice.
- Metric construction remains on demand and owns a further validated copy.

```mermaid
flowchart LR
    T[Path parameter in 0 to 1] --> P[UnitVolumeBainPath3D]
    P --> R[UnitVolumeBainPathPoint3D]
    R --> L[DirectLattice3D]
    R --> M[DirectLatticeMetric]
```

## Ownership boundary

Crystallographic mechanics and scientific use remain distinct.

- PhysKit owns the path equation, endpoint convention, primitive-basis orientation,
  lattice construction, and invariant evidence.
- KSD or another consumer owns path sampling, campaign selection, thresholds,
  electronic response, diagnostics, artifacts, and acceptance.
- The unit-volume lattice is dimensionless; material scale and units are external.
- The path is a represented crystallographic model, not evidence of a material phase
  transition or transformation mechanism.

## Deferred crystallography

Related owner work remains intentionally separate.

- The diamond-decorated fractional basis and FCC-to-Bain coordinate bridge are a later
  structure-level tranche.
- The area-preserving triangular-to-square rhombic path is a later 2D tranche.
- Generic homogeneous-deformation records require explicit active/passive and
  left/right-action conventions before introduction.
- Completion of this path authorizes none of those successors.

## Local mapping

Source, tests, and provenance remain explicit.

- Code: `src/python/projectkoios/physkit/periodic/lattice/deformation_paths.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/deformation_paths/`
- [Parent package](../index.md)
- [Extraction provenance](../../../../../../provenance/ksdft2effmass-crystallographic-paths.md)
