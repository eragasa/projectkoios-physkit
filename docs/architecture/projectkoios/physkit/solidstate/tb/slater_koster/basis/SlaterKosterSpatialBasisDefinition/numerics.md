# `SlaterKosterSpatialBasisDefinition` numerics

The class accepts a display identifier, nonempty unique center tuple, and nonempty unique coordinate tuple.

Its defining numerical relation is

$$
\mathcal B=(b_0,b_1,\ldots,b_{N_b-1}),\qquad \boldsymbol{\tau}_{a,:}=\tau_{c_a}.
$$

Every coordinate must reference a declared center. Equality and order are exact; dimensions and labels do not substitute for basis compatibility.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
