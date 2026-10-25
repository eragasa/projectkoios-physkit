# Package `projectkoios.physkit.linalg.matrix_norms`

## Responsibility

This package evaluates matrix norms for canonical PhysKit matrix quantities and restores
the exact input unit around unit-free Optimization results. Raw NumPy and SciPy values
use Optimization directly and are not accepted by this package.

```text
matrix_norms/
├── __init__.py
├── base.py
├── frobenius.py
├── maximum_absolute_entry.py
└── spectral.py
```

The package contains no behavioral module-level helpers, static utility methods,
classmethod helpers, Helper/Utils containers, factories, or compatibility adapters.
Shared behavior is owned by meaningful classes as instance methods.

## Accepted representations

- `MatrixQuantity`
- `ComplexMatrixQuantity`
- `SparseMatrixQuantity`
- `ComplexSparseMatrixQuantity`

The retained unit may be `PhysicalUnit` or `Unitless`. Every result preserves the exact
unit object. No new quantity-input union alias is defined.

## Supported exports

The package `__init__.py` and `__all__` expose:

- `BaseMatrixQuantityNorm`
- `MatrixNorm`
- exact Optimization `MatrixNormEnum`
- `FrobeniusMatrixNorm`
- `MaximumAbsoluteEntryMatrixNorm`
- `SpectralMatrixNorm`
- `SparseSpectralNormOptionsEnum`
- exact Optimization `ScipyCsrSpectralNormControls`

Optimization's `MatrixNormABC` and concrete backend implementations retain their owner
imports and are not re-exported through PhysKit.

## Numerical boundary

- Dense quantities delegate their exact immutable magnitude arrays.
- Sparse quantities cross their defined fresh canonical `to_csr()` boundary.
- Sparse Frobenius and maximum-entry norms never densify.
- Sparse spectral evaluation requires explicit `complete` or `iterative` selection.
- Only `complete` explicitly materializes CSR.
- `iterative` never densifies or falls back.

## Pages

- [Module `base`](base/index.md)
- [Module `frobenius`](frobenius/index.md)
- [Module `maximum_absolute_entry`](maximum_absolute_entry/index.md)
- [Module `spectral`](spectral/index.md)
- [Parent package](../index.md)
