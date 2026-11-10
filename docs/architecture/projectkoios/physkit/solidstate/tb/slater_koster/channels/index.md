# Module `projectkoios.physkit.solidstate.tb.slater_koster.channels`

This module owns the stable seven-channel enumeration and the nondimensional angular
matrix basis for homopolar $sp^3s^*$ two-center hopping.

## Current contract

- `HomopolarSp3sStarIntegralChannelEnum` names all seven channels.
- `HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER` fixes their observable order.
- `HomopolarSp3sStarTwoCenterChannelBasisConstructor` maps one normalized bond
  direction to seven real $5\times5$ coefficient matrices.
- `HomopolarSp3sStarTwoCenterChannelBasis` retains the direction and exact ordered
  matrix tuple as one correlated result.
- Lookup requires an exact enum value and a named argument.

The matrices contain angular coefficients only. They carry the nondimensional unit from
`UnitSystem.NONDIMENSIONAL` and contain no energy parameters.

## Class architecture

- [`HomopolarSp3sStarIntegralChannelEnum`](HomopolarSp3sStarIntegralChannelEnum/index.md)
- [`HomopolarSp3sStarTwoCenterChannelBasis`](HomopolarSp3sStarTwoCenterChannelBasis/index.md)
- [`HomopolarSp3sStarTwoCenterChannelBasisConstructor`](HomopolarSp3sStarTwoCenterChannelBasisConstructor/index.md)

## Exclusions

The module does not choose integral values, material meaning, shell meaning, fitting
policy, onsite terms, overlap, spin, or Bloch phases.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/channels.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/channels/`
- [Parent package](../index.md)
