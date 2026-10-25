# `MaximumAbsoluteEntryMatrixNorm`

## Mathematical identity

`MaximumAbsoluteEntryMatrixNorm` evaluates

\[
\lVert A\rVert_{\max} = \max_{i,j}\lvert A_{ij}\rvert.
\]

The definition applies to rectangular real and complex matrices. The class identifies
`MatrixNormEnum.maximum_absolute_entry`.

## Construction

The positional-only matrix argument accepts one canonical real/complex dense/sparse
matrix quantity. The inherited base retains the exact quantity object and rejects raw
numerical inputs.

## Dense execution

Dense quantities pass their exact immutable magnitude ndarray to
`NumpyMaximumAbsoluteEntryNorm(matrix).result()`. Complex entries use their absolute
magnitude rather than real-part ordering.

## Sparse execution

Sparse quantities pass one fresh canonical CSR value to
`ScipyCsrMaximumAbsoluteEntryNorm(matrix).result()`. The owner evaluates stored values
and the represented implicit-zero set without calling `toarray()`.

An all-zero sparse matrix therefore has norm `0.0`, including a canonical CSR value with
no stored entries.

## Result contract

The result is `ScalarQuantity(value, matrix.unit)` with exact unit identity. The scalar
magnitude is a built-in finite nonnegative `float`.

## Failure contract

Quantity construction owns finite represented input. Optimization failures propagate
unchanged. PhysKit performs no retry, fallback, coercion, or exception translation.
