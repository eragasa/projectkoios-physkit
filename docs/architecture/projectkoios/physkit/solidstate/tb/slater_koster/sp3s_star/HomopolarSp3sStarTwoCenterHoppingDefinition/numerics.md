# `HomopolarSp3sStarTwoCenterHoppingDefinition` numerics

The class accepts a `SlaterKosterBondDirection` and `HomopolarSp3sStarTwoCenterIntegrals`.

Its defining numerical relation is

$$
\mathcal D_T=(d,V).
$$

Both values require exact nominal owners. The definition performs no geometric or material inference.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
