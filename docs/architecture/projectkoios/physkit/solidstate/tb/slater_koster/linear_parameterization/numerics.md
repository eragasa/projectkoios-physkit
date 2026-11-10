# Linear Slater–Koster parameterization numerics

For parameter group $g$ and channel $\alpha$, each bond naming $g$ contributes the
nondimensional angular matrix $B_\alpha(d_b)$ to the appropriate source-target block.
The assembler also generates its transpose reverse. The resulting coefficient models
are

$$
B_{g,\alpha}(R)\in\mathbb C^{N_b\times N_b},
$$

stored as `complex128` with a nondimensional unit. For physical energy coefficients
$V_{g,\alpha}$ in one exact unit,

$$
H_R(V)=\sum_{g,\alpha}V_{g,\alpha}B_{g,\alpha}(R).
$$

The outer order is the declared parameter-group order. The inner order is

$$
(ss\sigma,sp\sigma,pp\sigma,pp\pi,s^*p\sigma,ss^*\sigma,s^*s^*\sigma).
$$

Equal displacements are accumulated before immutable terms are constructed and sorted
lexicographically. Nonfinite accumulation raises `OverflowError`; no clipping or
fallback occurs.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
