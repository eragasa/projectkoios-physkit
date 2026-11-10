# Class `SlaterKosterSpatialBasisCoordinate`

## Current contract

Identifies one ordered spin-independent spatial coordinate by center, angular component, and radial channel.

The public boundary accepts a center identifier, exact `RealCubicOrbitalComponentEnum`, and nonnegative built-in integer radial channel. Booleans, NumPy integers, negative channels, empty center labels, and arbitrary component strings are rejected.

The central represented relation is

$$
b_a=(c_a,\lambda_a,r_a).
$$

## Ownership boundary

The radial channel is coordinate identity and is not inferred to be a principal quantum number. In the current specialization, s star is `(s, 1)`. It does not represent spin, occupation, normalization, energy, or orbital wavefunctions.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/basis.py::SlaterKosterSpatialBasisCoordinate`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/basis/`
- [Owning module](../index.md)
