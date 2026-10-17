# Class `projectkoios.physkit.periodic.lattice.bain.path.UnitVolumeBainPathPoint3D`

## Current contract

The immutable record correlates one path parameter with its complete Bain geometry.

- Scalar fields retain the path parameter, axial ratio, and transverse/longitudinal
  scales.
- `direct_lattice` is an owned canonical `DirectLattice3D` copy.
- Construction verifies every scalar against the path law exactly.
- Construction verifies every direct-basis component against the declared orientation.
- `primitive_volume` delegates to the correlated direct lattice.

## Derived metric

The point exposes current metric and reciprocal geometry without duplicating state.

- `direct_metric` constructs `DirectLatticeMetric` on demand.
- The metric owns its own validated lattice copy, so object identity is not correlation
  evidence.
- Primitive-basis equality establishes point/metric correlation.
- Metric construction does not add physical units or material identity.

## Invariants

Invalid or ambiguous point state is rejected at construction.

- The path parameter is finite and lies in `[0, 1]`.
- Axial ratio and both scales are positive finite built-in floats.
- Scalar values and primitive vectors describe the same mathematical path point.
- Stored direct-basis arrays are immutable binary64 values.

## Local mapping

Source and tests remain adjacent to the module owner.

- Code: `src/python/projectkoios/physkit/periodic/lattice/bain/path.py::UnitVolumeBainPathPoint3D`
- Tests: `tests/projectkoios/physkit/periodic/lattice/bain/path/test__UnitVolumeBainPathPoint3D__contract.py`
- [Parent module](../index.md)
