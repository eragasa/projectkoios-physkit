# `MatrixNorm`

## Responsibility

`MatrixNorm` is the immutable, slotted public facade for semantic matrix-norm dispatch.
It validates PhysKit-facing identities and delegates to exactly one concrete
quantity-aware evaluator.

## Construction

```python
MatrixNorm(
    matrix,
    matrix_norm,
    *,
    sparse_spectral_option=None,
    sparse_spectral_controls=None,
)
```

- `matrix` is one canonical real/complex dense/sparse matrix quantity.
- `matrix_norm` is an exact Optimization `MatrixNormEnum` member.
- `sparse_spectral_option` is an exact `SparseSpectralNormOptionsEnum` member or `None`.
- `sparse_spectral_controls` is an exact `ScipyCsrSpectralNormControls` instance or
  `None`.

The sparse arguments are keyword-only so every sparse spectral policy choice remains
visible at the call site.

## Dispatch table

| `matrix_norm` | Quantity-aware implementation |
|---|---|
| `frobenius` | `FrobeniusMatrixNorm` |
| `spectral` | `SpectralMatrixNorm` |
| `maximum_absolute_entry` | `MaximumAbsoluteEntryMatrixNorm` |

The implementation map and diagnostic-name map are immutable internal state. No public
mapping accessor exists, and strings are not accepted as semantic selectors.

## Validation

Construction raises `TypeError` for a nonquantity matrix or a noncanonical enum,
option, or controls identity. Frobenius and maximum-entry construction raises
`ValueError` when supplied either sparse spectral argument. `SpectralMatrixNorm` owns
the representation-specific sparse argument decision table.

## Result lifecycle

`result()` constructs the selected evaluator and invokes its `.result()` exactly once.
The facade performs no numerical arithmetic, dense conversion, unit conversion, result
caching, retry, fallback, or exception translation.

## Immutability

The dataclass is frozen, slotted, and identity-based (`eq=False`). Public construction
state cannot be reassigned. Internal dispatch mappings use `MappingProxyType` and cannot
be mutated through the facade.
