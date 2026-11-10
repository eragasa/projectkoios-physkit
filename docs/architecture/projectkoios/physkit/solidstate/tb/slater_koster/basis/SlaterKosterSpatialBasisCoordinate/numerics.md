# `SlaterKosterSpatialBasisCoordinate` numerics

The class accepts a center identifier, exact `RealCubicOrbitalComponentEnum`, and nonnegative built-in integer radial channel.

Its defining numerical relation is

$$
b_a=(c_a,\lambda_a,r_a).
$$

Booleans, NumPy integers, negative channels, empty center labels, and arbitrary component strings are rejected.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
