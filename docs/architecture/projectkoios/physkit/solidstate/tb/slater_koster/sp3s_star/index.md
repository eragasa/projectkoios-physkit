# Module `projectkoios.physkit.solidstate.tb.slater_koster.sp3s_star`

This module owns seven physical two-center integrals and construction of one real
homopolar $sp^3s^*$ hopping matrix.

## Current contract

- `HomopolarSp3sStarTwoCenterIntegrals` retains seven real scalar energy quantities in
  stable channel order and one explicit physical `UnitSystem`.
- `HomopolarSp3sStarTwoCenterHoppingDefinition` binds those integrals to one normalized
  source-to-target bond direction.
- `HomopolarSp3sStarTwoCenterHoppingConstructor` forms the exact linear combination of
  the seven nondimensional angular matrices.
- `HomopolarSp3sStarTwoCenterHopping` correlates the definition, immutable matrix, unit,
  and fixed orbital order.
- Finite-input arithmetic that produces nonfinite entries raises `OverflowError`.

## Class architecture

- [`HomopolarSp3sStarTwoCenterIntegrals`](HomopolarSp3sStarTwoCenterIntegrals/index.md)
- [`HomopolarSp3sStarTwoCenterHoppingDefinition`](HomopolarSp3sStarTwoCenterHoppingDefinition/index.md)
- [`HomopolarSp3sStarTwoCenterHopping`](HomopolarSp3sStarTwoCenterHopping/index.md)
- [`HomopolarSp3sStarTwoCenterHoppingConstructor`](HomopolarSp3sStarTwoCenterHoppingConstructor/index.md)

## Exclusions

The module does not own lattice geometry, neighbor inventories, onsite terms, assembly,
Bloch phases, fitting, evidence, or scientific acceptance.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star/`
- [Parent package](../index.md)
