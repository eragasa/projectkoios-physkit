# Module `projectkoios.physkit.linalg.matrix.norm.frobenius`

## Responsibility

This module owns [`FrobeniusMatrixNorm`](FrobeniusMatrixNorm/index.md), the
quantity-aware adapter for the matrix Frobenius norm. It owns representation dispatch
and unit restoration and delegates all reduction arithmetic to Optimization.

## Numerical owners

- Dense real and complex quantities use
  `projectkoios.numerics.linalg.numpy.frobenius.NumpyFrobeniusNorm`.
- Sparse real and complex quantities use
  `projectkoios.numerics.linalg.scipy.frobenius.ScipyCsrFrobeniusNorm`.

The sparse path remains CSR. The module contains no dense fallback, local reduction,
Frobenius Power API, or squared-unit behavior.

## Related pages

- [Class `FrobeniusMatrixNorm`](FrobeniusMatrixNorm/index.md)
- [Matrix-norm package](../index.md)
- [Base contracts](../base/index.md)
