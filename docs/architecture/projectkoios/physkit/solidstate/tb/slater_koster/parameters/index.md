# Module `projectkoios.physkit.solidstate.tb.slater_koster.parameters`

This module owns explicit parameter-sharing groups for homopolar $sp^3s^*$ bonds.

## Current contract

`HomopolarSp3sStarParameterGroup` binds a nonempty identifier to one validated
`HomopolarSp3sStarTwoCenterIntegrals` value. The record is frozen, slotted, and
keyword-only.

A group identifier denotes only explicit coefficient sharing. Clients may attach a
material or shell interpretation, but that interpretation is not inferred or stored by
this reusable owner.

## Class architecture

- [`HomopolarSp3sStarParameterGroup`](HomopolarSp3sStarParameterGroup/index.md)

## Exclusions

The module does not define parameter databases, default materials, neighbor shells,
distance scaling, fitting, regularization, uncertainty, evidence, or acceptance.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/parameters.py`
- Tests: exercised through `tests/projectkoios/physkit/solidstate/tb/slater_koster/`
- [Parent package](../index.md)
