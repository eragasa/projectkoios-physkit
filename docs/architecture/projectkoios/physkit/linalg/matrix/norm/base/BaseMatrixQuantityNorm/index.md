# `BaseMatrixQuantityNorm`

## Responsibility

`BaseMatrixQuantityNorm` is the nominal abstract base for one matrix-quantity norm. It
owns representation-independent quantity behavior and delegates the unit-free magnitude
to one concrete mathematical norm.

## Construction

The positional-only `matrix` argument accepts exactly:

- `MatrixQuantity`;
- `ComplexMatrixQuantity`;
- `SparseMatrixQuantity`; or
- `ComplexSparseMatrixQuantity`.

Construction retains the exact object. A raw ndarray, SciPy sparse object, unrelated
quantity, or structurally similar object raises `TypeError`.

## Abstract contract

Concrete subclasses implement:

```python
@property
def matrix_norm(self) -> MatrixNormEnum: ...

def _evaluate_magnitude(self) -> float: ...
```

`matrix_norm` returns the exact Optimization semantic enum member.
`_evaluate_magnitude()` invokes one exact Optimization evaluator and returns its
built-in finite nonnegative result.

## Result lifecycle

`result()` performs one evaluation per call and constructs:

```python
ScalarQuantity(self._evaluate_magnitude(), self.matrix.unit)
```

The method does not cache results. It preserves exact quantity and unit identity and
performs no unit conversion.

## Representation access

- `_dense_matrix()` returns the exact immutable magnitude ndarray owned by a dense
  quantity and rejects a sparse quantity.
- `_csr_matrix()` calls `to_csr()` on a sparse quantity and returns the fresh canonical
  CSR value defined by the quantity contract.
- `_matrix_is_sparse()` performs nominal dense/CSR dispatch.

These internal methods centralize representation handling without introducing a raw
matrix facade or a replacement input alias.

## Ownership boundary

The class is distinct from Optimization's `MatrixNormABC`. Optimization's base returns a
unit-free `float`; the PhysKit base returns `ScalarQuantity`. Neither base subclasses
nor adapts the other.
