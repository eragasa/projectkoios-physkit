# `HomopolarSp3sStarTwoCenterHoppingConstructor` numerics

The class accepts one keyword-only `HomopolarSp3sStarTwoCenterHoppingDefinition`.

Its defining numerical relation is

$$
T(d;V)=\sum_{\alpha=1}^{7}V_\alpha B_\alpha(d),\qquad T(-d;V)=T(d;V)^{\mathsf T}.
$$

Binary64 accumulation is checked for nonfinite output. Overflow raises rather than clipping or silently regularizing values.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
