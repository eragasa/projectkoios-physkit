# `RealCubicOrbitalComponentEnum` numerics

The class accepts one exact enum member rather than an arbitrary string.

Its defining numerical relation is

$$
\lambda\in\{s,p_x,p_y,p_z,d_{xy},d_{yz},d_{zx},d_{x^2-y^2},d_{3z^2-r^2}\}.
$$

Enum membership supplies stable string values for exact records and serialization-facing labels.

All arrays remain inside the nominal PhysKit quantity owners documented by the module.
Tuple order, matrix shape, units, and enum order are observable contract data; the class
performs no hidden conversion, sorting, fallback, or string-based inference unless the
defining relation explicitly states otherwise.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
