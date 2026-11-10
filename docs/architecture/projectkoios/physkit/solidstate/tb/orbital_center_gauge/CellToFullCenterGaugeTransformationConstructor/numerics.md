# `CellToFullCenterGaugeTransformationConstructor` numerics

The class accepts one keyword-only `CellToFullCenterGaugeTransformationDefinition`.

Its defining numerical relation is

$$
D_{aa}(q)=e^{+2\pi i q\cdot\tau_a}.
$$

The action performs a rowwise center-wave-vector contraction, vectorized complex exponential, and diagonal embedding while preserving input order.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
