# Class `RealCubicOrbitalComponentEnum`

## Current contract

Defines the closed real-cubic angular-component vocabulary admitted by generic spatial-basis coordinates.

The public boundary accepts one exact enum member rather than an arbitrary string. Enum membership supplies stable string values for exact records and serialization-facing labels.

The central represented relation is

$$
\lambda\in\{s,p_x,p_y,p_z,d_{xy},d_{yz},d_{zx},d_{x^2-y^2},d_{3z^2-r^2}\}.
$$

## Ownership boundary

An admitted component name identifies a coordinate label; it does not assert that a hopping table or physical model exists for that component. It does not include spin, total-angular-momentum labels, complex spherical harmonics, or f orbitals.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/basis.py::RealCubicOrbitalComponentEnum`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/basis/`
- [Owning module](../index.md)
