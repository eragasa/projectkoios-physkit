# Module `projectkoios.physkit.solidstate.tb.slater_koster.linear_parameterization`

This module owns the reusable linear ordinary-cell representation of an admitted
homopolar $sp^3s^*$ bond inventory.

## Current contract

- `HomopolarSp3sStarLinearParameterizationAssembler` creates one nondimensional
  `DenseBlockHoppingModel` for every explicit parameter-group and stable channel pair.
- Output order is parameter-group major and channel minor.
- All seven channels remain present even when a concrete coefficient is zero.
- `HomopolarSp3sStarLinearParameterization` retains the exact source definition and
  supports exact group/channel lookup.
- Each model preserves basis identity, basis size, ordinary-cell displacement support,
  generated Hermitian reverses, and the source energy-reference identity.

## Class architecture

- [`HomopolarSp3sStarLinearParameterization`](HomopolarSp3sStarLinearParameterization/index.md)
- [`HomopolarSp3sStarLinearParameterizationAssembler`](HomopolarSp3sStarLinearParameterizationAssembler/index.md)

## Exclusions

The module does not own onsite terms, parameter fitting, design-matrix flattening,
residuals, sample selection, material semantics, evidence, or acceptance.

- [Numerics](numerics.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
- [Parent package](../index.md)
