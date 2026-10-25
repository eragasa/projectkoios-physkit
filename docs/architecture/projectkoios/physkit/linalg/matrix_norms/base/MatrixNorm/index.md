# `MatrixNorm`

`MatrixNorm` is the immutable quantity-aware facade for one selected matrix norm.

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

- `matrix` is a canonical real/complex dense/sparse matrix quantity.
- `matrix_norm` is the exact Optimization `MatrixNormEnum`.
- sparse spectral options are accepted only for sparse spectral evaluation.

`result()` returns one `ScalarQuantity` from the selected concrete quantity-aware class.

## Internal mappings

The facade owns immutable internal/non-public mappings from:

- norm enum to stable diagnostic name; and
- norm enum to quantity-aware implementation class.

Mappings use `MappingProxyType`, are not exported, and have no public accessor. Strings
are not selectors. `MatrixNormEnum.name` is the supported textual semantic identity.

## Dispatch

- `frobenius` selects `FrobeniusMatrixNorm`.
- `spectral` selects `SpectralMatrixNorm`.
- `maximum_absolute_entry` selects `MaximumAbsoluteEntryMatrixNorm`.

The facade performs no numerical calculation, dense materialization, unit conversion,
fallback, retry, caching, or exception translation.
