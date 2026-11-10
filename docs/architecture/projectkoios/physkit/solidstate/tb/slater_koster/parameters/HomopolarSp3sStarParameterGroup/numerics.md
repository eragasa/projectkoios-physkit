# `HomopolarSp3sStarParameterGroup` numerics

The class accepts a nonempty built-in string identifier and exact `HomopolarSp3sStarTwoCenterIntegrals` value.

Its defining numerical relation is

$$
g=(\operatorname{id}_g,V_g).
$$

The record validates nominal ownership only; uniqueness and use of identifiers require the complete model-definition inventory.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
