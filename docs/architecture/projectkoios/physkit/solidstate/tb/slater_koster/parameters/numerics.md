# Slater–Koster parameter-group numerics

A parameter group is the immutable pair

$$
g=(\operatorname{id}_g,V_g),
$$

where the integral vector is

$$
V_g=
\left(
V_{ss\sigma},
V_{sp\sigma},
V_{pp\sigma},
V_{pp\pi},
V_{s^*p\sigma},
V_{ss^*\sigma},
V_{s^*s^*\sigma}
\right)_g.
$$

Every component has one exact physical energy unit $U_E$:

$$
[V_{g,\alpha}]=U_E
\qquad\text{for all }\alpha.
$$

A parameterized bond carries the group identifier, and model-definition validation
resolves it by exact equality:

$$
g(b)=g
\quad\Longleftrightarrow\quad
\operatorname{id}_{g(b)}=\operatorname{id}_g.
$$

There is no numerical interpolation, conversion, fallback, or nearest-match lookup in
this module. Duplicate identifiers are rejected by the model-definition owner, where
the complete group inventory is available.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
