# Linear Slater–Koster Bloch numerics

For each nondimensional channel model,

$$
B_{g,\alpha}(q_j)=\sum_R B_{g,\alpha}(R)
\exp(+2\pi i q_j\cdot R).
$$

The positive exponent, integer ordinary-cell displacement, no-normalization convention,
and reciprocal-coordinate interpretation are inherited from
`DenseBlockHoppingBlochEvaluator`. Results retain nested order

$$
[g][\alpha][j].
$$

For energy coefficients $V_{g,\alpha}$,

$$
H_{\mathrm{cell}}(q_j;V)=
\sum_{g,\alpha}V_{g,\alpha}B_{g,\alpha}(q_j).
$$

Because only integer $R$ enters the phase, each channel matrix is elementwise periodic
under $q\mapsto q+G$ for integer reciprocal coordinate $G$. Orbital-center gauge
transformations remain a separate owner.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
