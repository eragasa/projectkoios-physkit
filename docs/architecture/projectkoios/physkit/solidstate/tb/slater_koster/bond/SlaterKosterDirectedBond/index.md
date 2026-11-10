# Class `SlaterKosterDirectedBond`

## Current contract

Identifies one directed three-dimensional ordinary-cell bond independently of physical Cartesian geometry.

The public boundary accepts one 3D `LatticeDisplacement` and nonempty source and target center identifiers. The lattice dimension must be three. A zero-displacement self block is rejected because onsite terms have separate ownership.

The central represented relation is

$$
b=(R,s,t),\qquad b^{\mathrm{rev}}=(-R,t,s).
$$

## Ownership boundary

The record controls ordinary-cell coefficient placement and Bloch phase only; it does not claim that R determines the physical bond vector. It does not retain bond length, direction, shell, species, or parameter values.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/bond.py::SlaterKosterDirectedBond`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/bond/`
- [Owning module](../index.md)
