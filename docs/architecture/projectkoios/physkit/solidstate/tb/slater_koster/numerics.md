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

Focused tests compare all channel matrices with independently assembled arrays,
exercise scales from `1e-150` through `1e150`, and verify
$T(-d)=T(d)^T$.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
