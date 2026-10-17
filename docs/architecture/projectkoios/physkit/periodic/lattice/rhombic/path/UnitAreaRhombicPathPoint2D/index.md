# Class `projectkoios.physkit.periodic.lattice.rhombic.path.UnitAreaRhombicPathPoint2D`

## Current contract

The immutable record correlates one path parameter with complete rhombic geometry.

- Scalar fields retain the path parameter, interior angle, and common basis scale.
- `direct_lattice` is an owned canonical `DirectLattice2D` copy.
- Construction verifies every scalar against the path law exactly.
- Construction verifies every direct-basis component against the declared orientation.
- `primitive_area` delegates to the correlated direct lattice.

## Derived metric

The point exposes metric and reciprocal geometry without duplicating state.

- `direct_metric` constructs `DirectLatticeMetric` on demand.
- The metric owns a validated lattice copy, so object identity is not correlation
  evidence.
- The covariant metric is `[[1, cos(theta)], [cos(theta), 1]] / sin(theta)`.
- Metric construction adds no physical scale or material identity.

## Invariants

Invalid or ambiguous point state is rejected at construction.

- The path parameter is finite and lies in `[0, 1]`.
- Angle and scale are positive finite built-in floats.
- Scalars and primitive vectors describe the same mathematical point.
- Stored direct-basis arrays are immutable binary64 values.

## Local mapping

Source and tests remain adjacent to the module owner.

- Code: `src/python/projectkoios/physkit/periodic/lattice/rhombic/path.py::UnitAreaRhombicPathPoint2D`
- Tests: `tests/projectkoios/physkit/periodic/lattice/rhombic/path/test__UnitAreaRhombicPathPoint2D__contract.py`
- [Parent module](../index.md)
