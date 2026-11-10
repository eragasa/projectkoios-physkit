# Class `HomopolarSp3sStarTwoCenterHoppingDefinition`

## Current contract

Binds one source-to-target bond direction to one exact seven-integral set.

The public boundary accepts a `SlaterKosterBondDirection` and `HomopolarSp3sStarTwoCenterIntegrals`. Both values require exact nominal owners. The definition performs no geometric or material inference.

The central represented relation is

$$
\mathcal D_T=(d,V).
$$

## Ownership boundary

The pair fully specifies the represented two-center matrix under the module convention, but not a lattice model or scientific fit. It does not retain cell displacement, center identifiers, onsite terms, or evidence.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star.py::HomopolarSp3sStarTwoCenterHoppingDefinition`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/sp3s_star/`
- [Owning module](../index.md)
