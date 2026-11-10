# `CellToFullCenterGaugeTransformation` numerics

The class accepts a validated transformation definition and one square `ComplexMatrixQuantity`.

Its defining numerical relation is

$$
D(q)=\operatorname{diag}\!\left(e^{+2\pi i q\cdot\tau_0},\ldots,e^{+2\pi i q\cdot\tau_{N_b-1}}\right).
$$

The matrix shape must equal `(basis_size, basis_size)` and its unit must be the nondimensional unit-system view.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
