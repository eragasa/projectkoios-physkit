# `HomopolarSp3sStarTwoCenterIntegrals` numerics

The class accepts seven `ScalarQuantity` values and one explicit non-nondimensional `UnitSystem`.

Its defining numerical relation is

$$
V=(V_{ss\sigma},V_{sp\sigma},V_{pp\sigma},V_{pp\pi},V_{s^*p\sigma},V_{ss^*\sigma},V_{s^*s^*\sigma}).
$$

Every scalar must carry a physical energy unit, all units must match exactly, and that unit must equal the selected system energy unit.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
