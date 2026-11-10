# `SlaterKosterDirectedBond` numerics

The class accepts one 3D `LatticeDisplacement` and nonempty source and target center identifiers.

Its defining numerical relation is

$$
b=(R,s,t),\qquad b^{\mathrm{rev}}=(-R,t,s).
$$

The lattice dimension must be three. A zero-displacement self block is rejected because onsite terms have separate ownership.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
