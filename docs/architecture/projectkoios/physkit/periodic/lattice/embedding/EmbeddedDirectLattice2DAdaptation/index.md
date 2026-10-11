# Class `EmbeddedDirectLattice2DAdaptation`

A frozen, slotted `ResultsObject` correlating:

- the exact request and original ambient `3×2` basis;
- canonical intrinsic `2×2` basis and `DirectLatticeMetric`;
- unit-aware induced metric and inverse metric;
- ambient dual reciprocal basis;
- unitless ordered orthonormal frame and tangent projector; and
- explicit fractional/ambient and intrinsic/ambient coordinate maps.

Public construction validates shapes, units, orthonormality, projector
idempotence, induced-metric agreement, duality, and intrinsic-metric agreement.
Methods named `project_ambient_...` explicitly project away a normal component.

- [Module index](../index.md)
- Code: `src/python/projectkoios/physkit/periodic/lattice/embedding/__init__.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/embedding/test__EmbeddedDirectLattice2DAdaptation__coordinate_maps.py`
