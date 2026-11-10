# Class `HomopolarSp3sStarParameterGroup`

## Current contract

Identifies one explicitly shared immutable seven-integral parameter set.

The public boundary accepts a nonempty built-in string identifier and exact `HomopolarSp3sStarTwoCenterIntegrals` value. The record validates nominal ownership only; uniqueness and use of identifiers require the complete model-definition inventory.

The central represented relation is

$$
g=(\operatorname{id}_g,V_g).
$$

## Ownership boundary

Group identity expresses declared coefficient sharing and carries no intrinsic material, shell, symmetry, or acceptance semantics. It does not own fitting, uncertainty, databases, provenance, or parameter conversion.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/parameters.py::HomopolarSp3sStarParameterGroup`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/`
- [Owning module](../index.md)
