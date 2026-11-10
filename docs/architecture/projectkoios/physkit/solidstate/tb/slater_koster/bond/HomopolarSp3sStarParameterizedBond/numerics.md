# `HomopolarSp3sStarParameterizedBond` numerics

The class accepts a validated homopolar sp3s-star geometry and a built-in string group identifier.

Its defining numerical relation is

$$
b\longmapsto g(b).
$$

The record accepts only the exact geometry owner and a nonempty built-in string. Group resolution occurs in model-definition assembly.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
