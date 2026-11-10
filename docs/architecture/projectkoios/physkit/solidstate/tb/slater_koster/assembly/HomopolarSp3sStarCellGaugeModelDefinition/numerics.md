# `HomopolarSp3sStarCellGaugeModelDefinition` numerics

The class accepts a keyword-only model identifier, admitted basis, group tuple, bond tuple, and energy-reference string.

Its defining numerical relation is

$$
\mathcal D_{\mathrm{cell}}=(\mathcal B,\{V_g\},\{b\},E_{\mathrm{ref}}).
$$

Group identifiers and directed bonds are unique, explicit reverses are absent, endpoints belong to the basis, units agree, and every group is referenced.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
