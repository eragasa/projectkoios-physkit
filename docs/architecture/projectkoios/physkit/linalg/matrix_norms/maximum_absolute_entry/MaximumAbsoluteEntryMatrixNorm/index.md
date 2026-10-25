# `MaximumAbsoluteEntryMatrixNorm`

`MaximumAbsoluteEntryMatrixNorm` derives from `BaseMatrixQuantityNorm` and identifies
`MatrixNormEnum.maximum_absolute_entry`.

## Evaluation

- Dense quantities delegate exact magnitude arrays to
  `NumpyMaximumAbsoluteEntryNorm(matrix).result()`.
- Sparse quantities delegate fresh canonical CSR values to
  `ScipyCsrMaximumAbsoluteEntryNorm(matrix).result()`.

Implicit sparse zeros participate through the Optimization owner. Sparse storage is
never materialized as dense.

The result preserves the exact matrix-entry unit. Backend validation and numerical
failures propagate unchanged.
