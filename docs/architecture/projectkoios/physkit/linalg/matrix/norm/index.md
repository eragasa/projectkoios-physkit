# Package `projectkoios.physkit.linalg.matrix.norm`

## Responsibility

This package evaluates Frobenius, spectral, and maximum-absolute-entry norms for
canonical PhysKit matrix quantities. PhysKit owns nominal quantity acceptance, exact
unit restoration, semantic dispatch, and explicit sparse spectral path selection.
`projectkoios-optimization` owns every unit-free numerical reduction and decomposition.

Raw NumPy arrays and SciPy sparse objects are not PhysKit norm inputs. Raw numerical
callers import Optimization evaluators directly.

## Source structure

```text
norm/
├── __init__.py
├── base.py
├── frobenius.py
├── maximum_absolute_entry.py
└── spectral.py
```

- [`base`](base/index.md) owns the nominal quantity contract and dispatch facade.
- [`frobenius`](frobenius/index.md) owns quantity-aware Frobenius evaluation.
- [`maximum_absolute_entry`](maximum_absolute_entry/index.md) owns quantity-aware
  maximum-entry evaluation.
- [`spectral`](spectral/index.md) owns dense and explicitly selected sparse spectral
  evaluation.

## Accepted representations

Every evaluator accepts exactly one nominal quantity instance:

- `MatrixQuantity` for dense real binary64 entries;
- `ComplexMatrixQuantity` for dense complex128 entries;
- `SparseMatrixQuantity` for canonical real binary64 CSR entries; or
- `ComplexSparseMatrixQuantity` for canonical complex128 CSR entries.

The package defines no replacement union alias. Quantity constructors own dimensionality,
dtype, canonical CSR structure, finiteness, and immutable storage. Norm evaluators retain
the exact quantity object.

## Result and unit contract

Every public `.result()` returns:

```python
ScalarQuantity(owner_result, matrix.unit)
```

The scalar magnitude is a built-in finite nonnegative `float` established by the
Optimization owner. The result retains the exact input unit object. A `PhysicalUnit`
remains physical and a `Unitless` instance remains unitless. Norm evaluation performs no
unit conversion or rescaling.

All three norms carry the matrix-entry unit. Frobenius Power and its squared-unit
contract remain outside this package.

## Representation dispatch

| Quantity representation | Frobenius | Maximum entry | Spectral |
|---|---|---|---|
| Dense real or complex | NumPy owner | NumPy owner | Complete NumPy SVD owner |
| CSR real or complex | SciPy CSR owner | SciPy CSR owner | Explicit `complete` or `iterative` path |

Sparse Frobenius, sparse maximum-entry, and iterative sparse spectral evaluation never
materialize the represented matrix. Complete sparse spectral evaluation contains the
single explicit CSR-to-dense boundary.

## Supported exports

`projectkoios.physkit.linalg.matrix.norm.__all__` exposes:

- `BaseMatrixQuantityNorm`
- `FrobeniusMatrixNorm`
- `MatrixNorm`
- exact Optimization `MatrixNormEnum`
- `MaximumAbsoluteEntryMatrixNorm`
- exact Optimization `ScipyCsrSpectralNormControls`
- `SparseSpectralNormOptionsEnum`
- `SpectralMatrixNorm`

The package does not re-export Optimization backend evaluator classes or
`MatrixNormABC`. It provides no factory, string selector, callable shortcut,
compatibility adapter, raw-input alias, `ncv` alias, or mapping accessor.

## Failure ownership

- Nominal identity violations raise `TypeError` at the PhysKit boundary.
- Missing, incompatible, or irrelevant sparse spectral arguments raise `ValueError`.
- Complete sparse materialization rejects empty matrices before allocation and rejects
  nonfinite materialized entries with `FloatingPointError`.
- Optimization retains validation and failure ownership for each numerical algorithm.
- Backend nonconvergence propagates without retry, exception translation, or fallback.

## Behavioral ownership

Production behavior resides in cohesive class instance methods. The package contains no
behavioral module-level functions, static utility methods, classmethod helpers,
Helper/Utils containers, or structural protocols.

## Related pages

- [Parent package `matrix`](../index.md)
- [Package `projectkoios.physkit.linalg`](../../index.md)
- [KSD prior-art provenance](../../../../../../provenance/ksdft2effmass-matrix-norms.md)
