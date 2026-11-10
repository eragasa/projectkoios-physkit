# `CellToFullCenterGaugeTransformationDefinition` numerics

The class accepts an immutable nondimensional center matrix of shape `(basis_size, 3)` and an immutable nondimensional wave vector of shape `(3,)`.

Its defining numerical relation is

$$
\theta_a=q\cdot\tau_a,\qquad a=0,\ldots,N_b-1.
$$

Validation is nominal, shape-exact, and unit-exact. It performs no center sorting, basis inference, dimensional projection, or reciprocal-axis permutation.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
