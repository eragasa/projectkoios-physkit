# Module `projectkoios.physkit.linalg.matrix_norms.frobenius`

The module owns [`FrobeniusMatrixNorm`](FrobeniusMatrixNorm/index.md), the quantity-aware
Frobenius norm. It contains no Frobenius Power API and no local numerical kernel.

Dense and CSR reductions delegate to their exact Optimization owners. Sparse storage is
never materialized as dense.

## Related pages

- [Matrix-norm package](../index.md)
- [Base contracts](../base/index.md)
