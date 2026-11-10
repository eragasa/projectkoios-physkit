# Slater–Koster assembly numerics

Let $\mathcal B$ be an exact ordered basis with $N_b$ coordinates. For center $c$, let
$E_c\in\{0,1\}^{N_b\times5}$ embed its ordered five-orbital block, so

$$
E_c^{\mathsf T}E_c=I_5.
$$

For a directed representative

$$
b=(R_b,s_b,t_b,d_b,g_b),
$$

assembly computes $T_b=T(d_b;V_{g_b})$ and accumulates

$$
H_{R_b}\mathrel{+}=E_{s_b}T_bE_{t_b}^{\mathsf T}.
$$

It also generates the Hermitian reverse

$$
H_{-R_b}\mathrel{+}=E_{t_b}T_b^{\mathsf T}E_{s_b}^{\mathsf T}.
$$

Because the admitted integrals are real,

$$
T_b^{\dagger}=T_b^{\mathsf T}.
$$

If $R_b=0$ and $s_b\neq t_b$, both contributions enter $H_0$, giving the paired
source-target blocks in one coefficient matrix. A zero-cell self block is rejected by
bond validation.

For every retained displacement $R$, the mutable accumulator starts as

$$
H_R=0_{N_b\times N_b}\in\mathbb C^{N_b\times N_b}.
$$

Binary64 two-center blocks are accumulated into `complex128` matrices. If finite
accumulation produces a nonfinite entry, assembly raises `OverflowError`; no clipping or
fallback occurs. Final terms are sorted lexicographically by

$$
(R_1,R_2,R_3).
$$

The resulting ordinary-cell model is evaluated elsewhere according to

$$
H_{\mathrm{cell}}(q)=\sum_R H_R e^{+2\pi i q\cdot R}.
$$

Assembly itself does not evaluate this Fourier sum.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
