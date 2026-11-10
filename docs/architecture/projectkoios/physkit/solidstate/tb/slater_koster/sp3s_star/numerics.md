# Homopolar $sp^3s^*$ hopping numerics

Use orbital order

$$
(s,p_x,p_y,p_z,s^*)
$$

and normalized direction $d=(l,m,n)^{\mathsf T}$. The integral vector is

$$
V=(V_{ss\sigma},V_{sp\sigma},V_{pp\sigma},V_{pp\pi},
V_{s^*p\sigma},V_{ss^*\sigma},V_{s^*s^*\sigma}).
$$

The represented matrix is

$$
T(d;V)=
\begin{pmatrix}
V_{ss\sigma}
&V_{sp\sigma}d^{\mathsf T}
&V_{ss^*\sigma}\\
-V_{sp\sigma}d
&V_{pp\pi}I_3+(V_{pp\sigma}-V_{pp\pi})dd^{\mathsf T}
&-V_{s^*p\sigma}d\\
V_{ss^*\sigma}
&V_{s^*p\sigma}d^{\mathsf T}
&V_{s^*s^*\sigma}
\end{pmatrix}.
$$

Equivalently,

$$
T(d;V)=\sum_{\alpha=1}^{7}V_\alpha B_\alpha(d).
$$

All integrals must use the selected physical energy unit:

$$
[V_\alpha]=U_E=\operatorname{UnitSystem.energy\_unit},
\qquad U_E\neq 1.
$$

Construction uses an initially zero binary64 matrix and ordered accumulation. With
floating-point evaluation denoted by $\operatorname{fl}$,

$$
\widehat T=
\operatorname{fl}\!\left(\sum_{\alpha=1}^{7}V_\alpha B_\alpha\right).
$$

No clipping or rescaling is applied. If any entry of $\widehat T$ is nonfinite, the
constructor raises `OverflowError`. Otherwise the output is an immutable
`MatrixQuantity` with exact energy unit $U_E$.

For real homopolar integrals,

$$
T(-d;V)=T(d;V)^{\mathsf T}.
$$

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
