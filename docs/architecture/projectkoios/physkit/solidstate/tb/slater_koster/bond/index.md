# Module `projectkoios.physkit.solidstate.tb.slater_koster.bond`

This module owns directed ordinary-cell bond definitions and their explicit physical
Cartesian geometry.

## Current contract

- `SlaterKosterBondDirection` retains a normalized nondimensional three-component
  Cartesian direction.
- `SlaterKosterDirectedBond` retains one three-dimensional integer ordinary-cell
  displacement and exact source and target center identifiers.
- `HomopolarSp3sStarBondGeometry` binds that directed bond to a nonzero physical
  Cartesian displacement in the selected `UnitSystem`.
- `HomopolarSp3sStarParameterizedBond` binds the geometry to one explicit parameter-
  group identifier.
- Reversal negates the ordinary-cell displacement and swaps source and target centers.
- Public records require named arguments.

## Class architecture

- [`SlaterKosterBondDirection`](SlaterKosterBondDirection/index.md)
- [`SlaterKosterDirectedBond`](SlaterKosterDirectedBond/index.md)
- [`HomopolarSp3sStarBondGeometry`](HomopolarSp3sStarBondGeometry/index.md)
- [`HomopolarSp3sStarParameterizedBond`](HomopolarSp3sStarParameterizedBond/index.md)

## Exclusions

The module does not infer physical geometry from integer displacement, center labels,
neighbor-shell names, files, or lattice metadata. It does not select parameter values,
assemble matrices, add onsite terms, or establish material validity.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/bond.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/bond/`
- [Parent package](../index.md)
