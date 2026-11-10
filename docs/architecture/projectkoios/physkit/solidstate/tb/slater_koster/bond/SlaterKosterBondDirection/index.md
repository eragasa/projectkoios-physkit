# Class `SlaterKosterBondDirection`

## Current contract

Retains one normalized nondimensional Cartesian source-to-target bond direction.

The public boundary accepts one three-component nondimensional `VectorQuantity`. The Euclidean norm is checked with absolute tolerance `128 * eps`; shape and nondimensional unit are exact.

The central represented relation is

$$
d=(l,m,n),\qquad \lVert d\rVert_2=1.
$$

## Ownership boundary

Its components enter the represented formulas as direction cosines. The record does not own the physical bond length. It does not infer direction from integer displacement, center labels, or lattice geometry.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/bond.py::SlaterKosterBondDirection`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/bond/test__SlaterKosterBondDirection__init.py`
- [Owning module](../index.md)
