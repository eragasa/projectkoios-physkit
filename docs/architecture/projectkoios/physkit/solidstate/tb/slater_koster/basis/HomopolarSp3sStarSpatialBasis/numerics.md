# `HomopolarSp3sStarSpatialBasis` numerics

The class accepts one keyword-only `SlaterKosterSpatialBasisDefinition`.

Its defining numerical relation is

$$
\mathcal B_c=((c,s,0),(c,p_x,0),(c,p_y,0),(c,p_z,0),(c,s,1)),\qquad N_b=5N_c.
$$

The complete coordinate tuple must equal the generated expected tuple. Reordering, omission, insertion, or component substitution fails closed.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
