# `HomopolarSp3sStarBondGeometry` numerics

The class accepts a `SlaterKosterDirectedBond`, physical length `VectorQuantity`, and non-nondimensional `UnitSystem`.

Its defining numerical relation is

$$
d=\frac{\Delta r}{\lVert\Delta r\rVert_2},\qquad [\Delta r]=U_{\mathrm{system}}^{\mathrm{length}}.
$$

Length dimensionality, exact selected length unit, shape `(3,)`, and nonzero norm are validated before direction derivation.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
