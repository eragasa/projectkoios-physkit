# Module `projectkoios.physkit.linalg.matrix.norm.spectral`

## Responsibility

This module owns the quantity-aware spectral norm and the explicit selector for sparse
execution:

- [`SpectralMatrixNorm`](SpectralMatrixNorm/index.md)
- [`SparseSpectralNormOptionsEnum`](SparseSpectralNormOptionsEnum/index.md)

It owns dense-versus-sparse path validation, the explicit complete sparse materialization
boundary, iterative-control forwarding, and exact unit restoration. Optimization owns
complete SVD, direct small-shape sparse arithmetic, ARPACK iteration, residual
qualification, and convergence failure classification.

## Mathematical identity

The spectral norm is the largest singular value of a general rectangular matrix:

\[
\lVert A\rVert_2 = \sigma_{\max}(A).
\]

The module does not replace singular values with eigenvalue magnitudes.

## Numerical owners

- Dense quantities use `NumpySpectralNorm`.
- Sparse `complete` evaluation materializes CSR once and then uses
  `NumpySpectralNorm`.
- Sparse `iterative` evaluation uses `ScipyCsrSpectralNorm` with exact
  `ScipyCsrSpectralNormControls`.

## Related pages

- [Class `SpectralMatrixNorm`](SpectralMatrixNorm/index.md)
- [Enum `SparseSpectralNormOptionsEnum`](SparseSpectralNormOptionsEnum/index.md)
- [Matrix-norm package](../index.md)
- [Base contracts](../base/index.md)
