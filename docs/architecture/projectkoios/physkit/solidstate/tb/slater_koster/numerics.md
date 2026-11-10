# Slater–Koster numerics

For a unit direction $d=(l,m,n)$, the fixed orbital order is
$(s,p_x,p_y,p_z,s^*)$. The seven energy channels are ordered as

1. $V_{ss\sigma}$,
2. $V_{sp\sigma}$,
3. $V_{pp\sigma}$,
4. $V_{pp\pi}$,
5. $V_{s^*p\sigma}$,
6. $V_{ss^*\sigma}$,
7. $V_{s^*s^*\sigma}$.

The represented matrix is

$$
T(d;V)=\sum_\alpha V_\alpha B_\alpha(d),
$$

where every $B_\alpha$ is a nondimensional real $5\times5$ matrix. The nonzero blocks
are

$$
T_{ss}=V_{ss\sigma},\qquad
T_{s,p_i}=d_iV_{sp\sigma},\qquad
T_{p_i,s}=-d_iV_{sp\sigma},
$$

$$
T_{p_i,p_j}=d_i d_jV_{pp\sigma}
 +(\delta_{ij}-d_i d_j)V_{pp\pi},
$$

$$
T_{s^*,p_i}=d_iV_{s^*p\sigma},\qquad
T_{p_i,s^*}=-d_iV_{s^*p\sigma},
$$

$$
T_{s,s^*}=T_{s^*,s}=V_{ss^*\sigma},\qquad
T_{s^*,s^*}=V_{s^*s^*\sigma}.
$$

All seven coefficients use the energy unit selected by the owning `UnitSystem`.
The API value is a normalized Cartesian bond direction; its three components enter
these formulas as direction cosines. It has Euclidean norm one to a tolerance of
`128 * eps`. Construction uses binary64 arithmetic and raises `OverflowError` when
finite inputs produce nonfinite matrix entries. Values are not clipped or silently
regularized.

For a directed representative bond $b=(R_b,s_b,t_b,d_b,g_b)$, ordinary-cell
assembly embeds the five-orbital block as

$$
H_{R_b}\mathrel{+}=E_{s_b}T(d_b;V_{g_b})E_{t_b}^{\mathsf T}
$$

and generates the Hermitian reverse

$$
H_{-R_b}\mathrel{+}=E_{t_b}T(d_b;V_{g_b})^{\mathsf T}E_{s_b}^{\mathsf T}.
$$

For zero displacement between distinct centers, both contributions accumulate in one
coefficient matrix. Explicit reverse representatives are rejected. Coefficients are
stored in lexicographic displacement order with `complex128` arithmetic. Assembly
raises `OverflowError` if finite block accumulation becomes nonfinite.

Focused tests compare all channel matrices with independently assembled arrays,
exercise scales from `1e-150` through `1e150`, verify $T(-d)=T(d)^T$, and check exact
forward/reverse block placement in zero and nonzero cells.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
