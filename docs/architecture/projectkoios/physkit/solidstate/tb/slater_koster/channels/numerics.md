# Slater–Koster channel-basis numerics

Let

$$
d=(l,m,n)^{\mathsf T},
\qquad \lVert d\rVert_2=1,
$$

and use orbital order $(s,p_x,p_y,p_z,s^*)$. The channel matrices
$B_\alpha(d)\in\mathbb R^{5\times5}$ are defined so that

$$
T(d;V)=\sum_{\alpha=1}^{7}V_\alpha B_\alpha(d).
$$

With $I_3$ the $3\times3$ identity, the seven matrices are

$$
B_{ss\sigma}=
\begin{pmatrix}
1&0&0\\0&0_{3\times3}&0\\0&0&0
\end{pmatrix},
$$

$$
B_{sp\sigma}=
\begin{pmatrix}
0&d^{\mathsf T}&0\\-d&0_{3\times3}&0\\0&0&0
\end{pmatrix},
$$

$$
B_{pp\sigma}=
\begin{pmatrix}
0&0&0\\0&dd^{\mathsf T}&0\\0&0&0
\end{pmatrix},
\qquad
B_{pp\pi}=
\begin{pmatrix}
0&0&0\\0&I_3-dd^{\mathsf T}&0\\0&0&0
\end{pmatrix},
$$

$$
B_{s^*p\sigma}=
\begin{pmatrix}
0&0&0\\0&0_{3\times3}&-d\\0&d^{\mathsf T}&0
\end{pmatrix},
$$

$$
B_{ss^*\sigma}=
\begin{pmatrix}
0&0&1\\0&0_{3\times3}&0\\1&0&0
\end{pmatrix},
\qquad
B_{s^*s^*\sigma}=
\begin{pmatrix}
0&0&0\\0&0_{3\times3}&0\\0&0&1
\end{pmatrix}.
$$

Construction uses binary64 arrays. The matrix tuple is ordered exactly as

$$
(ss\sigma,sp\sigma,pp\sigma,pp\pi,s^*p\sigma,ss^*\sigma,s^*s^*\sigma).
$$

Each matrix is wrapped in an immutable nondimensional `MatrixQuantity`. Lookup is an
index operation in this fixed enum order.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
