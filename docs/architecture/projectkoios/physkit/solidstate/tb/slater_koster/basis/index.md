# Module `projectkoios.physkit.solidstate.tb.slater_koster.basis`

This module owns exact, ordered, spin-independent spatial-basis definitions for reusable
Slater–Koster representations.

## Current contract

- `SlaterKosterBasisCenter` binds a nonempty center identifier to one nondimensional
  primitive fractional coordinate.
- `SlaterKosterSpatialBasisCoordinate` binds a center, real-cubic orbital component,
  and nonnegative radial channel.
- `SlaterKosterSpatialBasisDefinition` retains the ordered center and coordinate
  inventories without inferring compatibility from labels or dimensions.
- `HomopolarSp3sStarSpatialBasis` admits only exact five-coordinate blocks ordered as
  `(s[0], p_x[0], p_y[0], p_z[0], s[1])` for every declared center.
- Public records and lookups require named arguments.

The conventional $s^*$ coordinate is represented as radial channel one of the $s$
component. It is not introduced as a separate universal angular component.

## Class architecture

- [`RealCubicOrbitalComponentEnum`](RealCubicOrbitalComponentEnum/index.md)
- [`SlaterKosterBasisCenter`](SlaterKosterBasisCenter/index.md)
- [`SlaterKosterSpatialBasisCoordinate`](SlaterKosterSpatialBasisCoordinate/index.md)
- [`SlaterKosterSpatialBasisDefinition`](SlaterKosterSpatialBasisDefinition/index.md)
- [`HomopolarSp3sStarSpatialBasis`](HomopolarSp3sStarSpatialBasis/index.md)

## Exclusions

The module does not select chemical species, lattice sites, neighbor shells, spin,
overlap, onsite energies, hopping coefficients, or scientific acceptance.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/basis.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/basis/`
- [Parent package](../index.md)
