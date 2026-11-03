# Block-hopping Bloch numerics

For fractional reciprocal coordinates $\mathbf q$ and integer displacements
$\mathbf R$, both routes evaluate

$$
H(\mathbf q)=\sum_{\mathbf R}T(\mathbf R)
\exp(+2\pi i\,\mathbf q\cdot\mathbf R).
$$

The exponent is positive and no normalization factor is applied. Dense evaluation
uses either exact ordered entrywise finite-series evaluation or the configured
shared-grid batched backend. Sparse evaluation accumulates phase-weighted CSR blocks,
canonicalizes the CSR result, and rejects nonfinite arithmetic. Results retain the
input coordinate record and model energy unit.

No route diagonalizes the matrix, checks Hermiticity, or applies a tolerance.

## Navigation

- [Module index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
