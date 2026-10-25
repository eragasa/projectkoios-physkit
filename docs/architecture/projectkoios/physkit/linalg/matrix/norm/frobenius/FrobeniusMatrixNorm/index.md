# `FrobeniusMatrixNorm`

## Mathematical identity

`FrobeniusMatrixNorm` evaluates

\[
\lVert A\rVert_F = \sqrt{\sum_{i,j}\lvert A_{ij}\rvert^2}.
\]

The definition applies to rectangular real and complex matrices. The class identifies
`MatrixNormEnum.frobenius`.

## Construction

The positional-only matrix argument accepts one canonical real/complex dense/sparse
matrix quantity. The inherited base retains the exact quantity object and rejects raw
numerical inputs.

## Dense execution

For `MatrixQuantity` and `ComplexMatrixQuantity`, the class passes the exact immutable
`matrix.magnitude` ndarray to:

```python
NumpyFrobeniusNorm(matrix.magnitude).result()
```

Optimization owns stable scaling, reduction order, representation checks, and finite
result enforcement.

## Sparse execution

For `SparseMatrixQuantity` and `ComplexSparseMatrixQuantity`, the class obtains one fresh
canonical CSR value through `matrix.to_csr()` and invokes:

```python
ScipyCsrFrobeniusNorm(csr).result()
```

The sparse path does not call `toarray()` and does not allocate a dense matrix.

## Result contract

The result is `ScalarQuantity(value, matrix.unit)`. The unit is the exact matrix-entry
unit object. The unit is not squared because the square root completes the norm.
Frobenius Power remains a separate Optimization capability.

## Failure contract

Quantity construction owns finite represented input. Optimization failures propagate
unchanged. PhysKit performs no retry, fallback, coercion, or exception translation.
