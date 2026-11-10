# `HomopolarSp3sStarTwoCenterChannelBasis` numerics

The class accepts a `SlaterKosterBondDirection` and a seven-element tuple of `MatrixQuantity` values.

Its defining numerical relation is

$$
\mathcal B_{\mathrm{channel}}(d)=(B_1(d),\ldots,B_7(d)),\qquad B_\alpha\in\mathbb R^{5\times5}.
$$

Every matrix has shape `(5, 5)` and the exact nondimensional unit. Lookup uses the stable enum order.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
