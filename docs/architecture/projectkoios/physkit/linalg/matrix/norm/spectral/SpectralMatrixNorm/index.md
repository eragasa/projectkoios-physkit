# `SpectralMatrixNorm`

## Mathematical identity

`SpectralMatrixNorm` evaluates the largest singular value

\[
\lVert A\rVert_2 = \sigma_{\max}(A)
\]

for rectangular real and complex matrix quantities. The class identifies
`MatrixNormEnum.spectral`.

## Construction

```python
SpectralMatrixNorm(
    matrix,
    *,
    sparse_spectral_option=None,
    sparse_spectral_controls=None,
)
```

The matrix is positional-only. Sparse policy is keyword-only and therefore explicit at
every call site.

## Argument decision table

| Representation | Option | Controls | Result |
|---|---|---|---|
| Dense | `None` | `None` | Complete dense evaluation |
| Dense | supplied | any | `ValueError` |
| Dense | `None` | supplied | `ValueError` |
| Sparse | `None` | any | `ValueError` |
| Sparse | `complete` | `None` | Explicit materialization and complete SVD |
| Sparse | `complete` | supplied | `ValueError` |
| Sparse | `iterative` | supplied | CSR owner evaluation |
| Sparse | `iterative` | `None` | `ValueError` |

A noncanonical option or controls object raises `TypeError` before this decision table is
applied.

## Dense execution

Dense quantities pass their exact immutable magnitude ndarray to
`NumpySpectralNorm(matrix).result()`. The Optimization owner validates the dense
representation and performs complete singular-value decomposition.

## Complete sparse execution

The `complete` path performs these operations in order:

1. obtain one fresh canonical CSR value through `matrix.to_csr()`;
2. reject a zero-length matrix axis with `ValueError`;
3. call `toarray()` exactly once under a local floating-warning boundary;
4. reject nonfinite materialized entries with `FloatingPointError`;
5. pass the resulting ndarray to `NumpySpectralNorm`; and
6. restore the exact matrix unit around the owner result.

This path is the package's only sparse dense-materialization boundary.

## Iterative sparse execution

The `iterative` path passes one fresh canonical CSR value and the exact retained controls
to:

```python
ScipyCsrSpectralNorm(csr, controls).result()
```

The Optimization owner uses a deterministic nondensifying small-Gram calculation for
`min(shape) <= 2` and ARPACK for `min(shape) >= 3`. Both branches belong to the same
selected sparse owner. PhysKit does not inspect shape to select a different algorithm.

The exact controls fields are:

- `tolerance`;
- `maximum_iterations`;
- `krylov_dimension`; and
- `seed`.

## Result and failure contract

The result is `ScalarQuantity(value, matrix.unit)` with exact unit identity. Dense SVD
and sparse iteration nonconvergence propagate as `RuntimeError` from Optimization.
Nonfinite owner arithmetic propagates as `FloatingPointError`. The class catches no
owner failure and performs no retry or fallback.
