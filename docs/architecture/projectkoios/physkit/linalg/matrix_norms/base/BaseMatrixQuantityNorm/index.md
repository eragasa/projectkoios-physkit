# `BaseMatrixQuantityNorm`

`BaseMatrixQuantityNorm` is the nominal PhysKit base for one quantity-aware matrix norm.

## Construction

Construction accepts exactly one of:

- `MatrixQuantity`
- `ComplexMatrixQuantity`
- `SparseMatrixQuantity`
- `ComplexSparseMatrixQuantity`

The exact quantity object is retained. Raw ndarray and CSR values are rejected because
raw numerical use belongs directly to Optimization.

## Contract

Concrete subclasses provide the canonical `MatrixNormEnum` identity and unit-free
magnitude evaluation. `result()` wraps the exact built-in finite nonnegative owner
result as:

```python
ScalarQuantity(value, matrix.unit)
```

The input unit object is retained whether it is physical or unitless.

The base owns shared instance behavior for dense magnitude access, canonical sparse
`to_csr()` access, representation classification, and unit restoration. It performs no
unit conversion or numerical reduction.

## Nominal boundary

This PhysKit base is distinct from Optimization's exact unit-free `MatrixNormABC`.
PhysKit does not redefine or subclass the Optimization ABC because its scalar result
contract returns `float`, not `ScalarQuantity`.
