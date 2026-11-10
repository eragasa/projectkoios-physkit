# `HomopolarSp3sStarTwoCenterHopping` numerics

The class accepts a hopping definition and `(5, 5)` energy-valued `MatrixQuantity`.

Its defining numerical relation is

$$
T(d;V)=\sum_{\alpha=1}^{7}V_\alpha B_\alpha(d).
$$

The matrix shape and energy unit must match the definition integrals exactly. The fixed orbital order is exposed explicitly.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
