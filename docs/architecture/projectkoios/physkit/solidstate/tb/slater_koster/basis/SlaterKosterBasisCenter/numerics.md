# `SlaterKosterBasisCenter` numerics

The class accepts a built-in string identifier and a nondimensional three-component `VectorQuantity`.

Its defining numerical relation is

$$
c\longmapsto\tau_c\in\mathbb R^3.
$$

The coordinate shape is exactly `(3,)`; no Cartesian lattice basis or unit conversion is accepted at this boundary.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
