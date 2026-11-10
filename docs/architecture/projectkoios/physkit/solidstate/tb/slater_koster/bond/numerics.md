# Slater–Koster bond numerics

A directed bond is the exact record

$$
b=(R,s,t),
\qquad R\in\mathbb Z^3,
$$

where $s$ and $t$ are explicit source and target center identifiers. Its Hermitian-
reverse directed bond is

$$
b^{\mathrm{rev}}=(-R,t,s).
$$

The bound physical geometry adds a nonzero Cartesian displacement

$$
\Delta r\in\mathbb R^3,
\qquad [\Delta r]=\text{length},
\qquad \lVert\Delta r\rVert_2>0.
$$

The normalized bond direction is derived only from that physical vector:

$$
d=(l,m,n)=\frac{\Delta r}{\lVert\Delta r\rVert_2},
\qquad \lVert d\rVert_2=1.
$$

`SlaterKosterBondDirection` checks the norm with absolute tolerance

$$
\varepsilon_d=128\,\epsilon_{64},
$$

where $\epsilon_{64}$ is binary64 machine epsilon. The direction is stored as a
`VectorQuantity` using `UnitSystem.NONDIMENSIONAL.length_unit`.

The physical displacement must use exactly the selected unit system's length unit:

$$
[\Delta r]=U_{\mathrm{system}}^{\mathrm{length}}.
$$

Unit compatibility alone is not enough at this boundary; conversion is an explicit
upstream action. A zero-cell self bond,

$$
R=0\quad\land\quad s=t,
$$

is rejected because onsite ownership is separate.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
