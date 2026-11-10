# Module `projectkoios.physkit.solidstate.tb.slater_koster.assembly`

This module owns ordinary-cell embedding of explicit homopolar $sp^3s^*$ directed-bond
representatives into PhysKit dense block-hopping models.

## Current contract

- `HomopolarSp3sStarCellGaugeModelDefinition` correlates one exact structured basis,
  unique parameter groups, unique directed representatives, and an explicit energy
  reference.
- Every bond endpoint must occur in the admitted basis, every bond must resolve one
  declared group, and every declared group must be used.
- Explicit Hermitian-reverse representatives are rejected because assembly generates
  them.
- `HomopolarSp3sStarCellGaugeModelAssembler` accumulates five-orbital blocks by exact
  integer displacement and emits lexicographically sorted `DenseBlockHoppingTerm`
  values.
- `HomopolarSp3sStarCellGaugeModelAssembly` retains the exact definition and generated
  `DenseBlockHoppingModel` as a correlated result.

## Class architecture

- [`HomopolarSp3sStarCellGaugeModelDefinition`](HomopolarSp3sStarCellGaugeModelDefinition/index.md)
- [`HomopolarSp3sStarCellGaugeModelAssembly`](HomopolarSp3sStarCellGaugeModelAssembly/index.md)
- [`HomopolarSp3sStarCellGaugeModelAssembler`](HomopolarSp3sStarCellGaugeModelAssembler/index.md)

## Exclusions

The module does not add onsite terms, generate lattices or neighbors, evaluate Bloch
matrices, construct linear fitting designs, fit parameters, retain evidence, or decide
scientific adequacy.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/assembly.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/assembly/`
- [Parent package](../index.md)
