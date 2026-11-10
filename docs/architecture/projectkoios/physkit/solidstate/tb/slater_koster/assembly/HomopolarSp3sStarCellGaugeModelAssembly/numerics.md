# `HomopolarSp3sStarCellGaugeModelAssembly` numerics

The class accepts a cell-gauge model definition and exact `DenseBlockHoppingModel`.

Its defining numerical relation is

$$
\mathcal A=(\mathcal D_{\mathrm{cell}},\{H_R\}_{R\in\mathbb Z^3}).
$$

Model identity, dimension, energy unit, energy reference, basis metadata, and basis size must agree with the retained definition.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
