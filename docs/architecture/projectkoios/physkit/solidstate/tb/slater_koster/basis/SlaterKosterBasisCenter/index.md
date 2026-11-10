# Class `SlaterKosterBasisCenter`

## Current contract

Binds one nonempty center identifier to an immutable primitive fractional coordinate.

The public boundary accepts a built-in string identifier and a nondimensional three-component `VectorQuantity`. The coordinate shape is exactly `(3,)`; no Cartesian lattice basis or unit conversion is accepted at this boundary.

The central represented relation is

$$
c\longmapsto\tau_c\in\mathbb R^3.
$$

## Ownership boundary

The center is an embedding label for represented orbitals, not an atom, species, Wannier center, or symmetry site unless a client establishes that meaning. It does not own lattice vectors, Cartesian positions, species, occupancy, or symmetry data.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/basis.py::SlaterKosterBasisCenter`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/basis/`
- [Owning module](../index.md)
