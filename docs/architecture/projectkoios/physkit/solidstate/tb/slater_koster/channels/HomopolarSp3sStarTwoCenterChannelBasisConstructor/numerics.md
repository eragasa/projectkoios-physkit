# `HomopolarSp3sStarTwoCenterChannelBasisConstructor` numerics

The class accepts one keyword-only `SlaterKosterBondDirection`.

Its defining numerical relation is

$$
B_\alpha(d)=\frac{\partial T(d;V)}{\partial V_\alpha}.
$$

The action uses exact channel dispatch and binary64 arrays; each branch writes only the entries owned by that channel.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
