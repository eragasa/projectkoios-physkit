# `FrobeniusMatrixNorm`

`FrobeniusMatrixNorm` derives from `BaseMatrixQuantityNorm` and identifies
`MatrixNormEnum.frobenius`.

## Evaluation

- `MatrixQuantity` and `ComplexMatrixQuantity` pass their exact magnitude arrays to
  `NumpyFrobeniusNorm(matrix).result()`.
- `SparseMatrixQuantity` and `ComplexSparseMatrixQuantity` pass a fresh canonical CSR
  value to `ScipyCsrFrobeniusNorm(matrix).result()`.

The class returns `ScalarQuantity(value, quantity.unit)` with exact unit identity. A norm
has the matrix-entry unit; the class does not expose Frobenius Power or square the unit.

## Boundary

Optimization owns representation validation, duplicate consolidation, arithmetic, and
finite-result enforcement. PhysKit does not catch or translate owner failures.
