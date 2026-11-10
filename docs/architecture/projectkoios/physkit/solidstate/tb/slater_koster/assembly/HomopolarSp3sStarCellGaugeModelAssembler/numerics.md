# `HomopolarSp3sStarCellGaugeModelAssembler` numerics

The class accepts one keyword-only `HomopolarSp3sStarCellGaugeModelDefinition`.

Its defining numerical relation is

$$
H_{R_b}\mathrel{+}=E_{s_b}T_bE_{t_b}^{\mathsf T},\qquad H_{-R_b}\mathrel{+}=E_{t_b}T_b^{\mathsf T}E_{s_b}^{\mathsf T}.
$$

The action resolves groups exactly, accumulates into complex128 matrices, raises on nonfinite accumulation, and sorts displacement terms lexicographically.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
