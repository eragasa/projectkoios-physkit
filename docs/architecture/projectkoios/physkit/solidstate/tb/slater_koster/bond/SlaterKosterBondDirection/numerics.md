# `SlaterKosterBondDirection` numerics

The class accepts one three-component nondimensional `VectorQuantity`.

Its defining numerical relation is

$$
d=(l,m,n),\qquad \lVert d\rVert_2=1.
$$

The Euclidean norm is checked with absolute tolerance `128 * eps`; shape and nondimensional unit are exact.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
