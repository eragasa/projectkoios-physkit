# `HomopolarSp3sStarIntegralChannelEnum` numerics

The class accepts one exact enum member rather than an arbitrary channel string.

Its defining numerical relation is

$$
\alpha\in(ss\sigma,sp\sigma,pp\sigma,pp\pi,s^*p\sigma,ss^*\sigma,s^*s^*\sigma).
$$

Enum values supply stable slugs used consistently by integral names and channel-matrix lookup.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
