# `SpectralMatrixNorm`

`SpectralMatrixNorm` derives from `BaseMatrixQuantityNorm` and identifies
`MatrixNormEnum.spectral`, the largest singular value of a general rectangular matrix.
It never substitutes an eigenvalue magnitude.

## Dense quantities

Dense quantities require sparse option and controls to be absent. Their exact magnitude
array is passed to `NumpySpectralNorm(matrix).result()`, which performs complete dense
SVD.

## Sparse quantities

Sparse quantities require one explicit `SparseSpectralNormOptionsEnum` member.

### Complete

The `complete` path:

1. forbids iterative controls;
2. obtains a fresh canonical CSR value from the quantity;
3. explicitly materializes one dense ndarray;
4. rejects nonfinite materialization arithmetic immediately; and
5. delegates to `NumpySpectralNorm`.

This is the only sparse densification path.

### Iterative

The `iterative` path:

1. requires an exact `ScipyCsrSpectralNormControls` instance;
2. obtains a fresh canonical CSR value;
3. delegates to `ScipyCsrSpectralNorm(matrix, controls).result()`; and
4. never densifies, retries, or falls back.

The Optimization evaluator uses a deterministic nondensifying direct small-Gram branch
for `min(shape) <= 2` and approximate ARPACK iteration for `min(shape) >= 3`. The direct
small-shape branch remains part of the selected sparse owner evaluator rather than a
third PhysKit option.

## Controls and failures

The public control field is `krylov_dimension`; `ncv` exists only at Optimization's
internal SciPy boundary. Owner validation, nonconvergence, nonfinite arithmetic, and
residual failures propagate without PhysKit translation.

The result preserves the exact matrix-entry unit whether physical or unitless.
