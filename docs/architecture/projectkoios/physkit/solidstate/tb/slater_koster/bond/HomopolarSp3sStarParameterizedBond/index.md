# Class `HomopolarSp3sStarParameterizedBond`

## Current contract

Binds one explicit bond geometry to one nonempty parameter-group identifier.

The public boundary accepts a validated homopolar sp3s-star geometry and a built-in string group identifier. The record accepts only the exact geometry owner and a nonempty built-in string. Group resolution occurs in model-definition assembly.

The central represented relation is

$$
b\longmapsto g(b).
$$

## Ownership boundary

The identifier states coefficient sharing only. It does not intrinsically encode material, shell, symmetry, or fit status. It does not contain integral values or resolve group existence by itself.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/bond.py::HomopolarSp3sStarParameterizedBond`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/bond/`
- [Owning module](../index.md)
