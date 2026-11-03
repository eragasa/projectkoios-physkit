# `DenseBlockHoppingTerm` numerics

The coefficient must be an exact immutable `ComplexMatrixQuantity`. Its shape is
reported through `coefficient_shape`; the abstract owner validates that it is
nonempty and square and derives `basis_size`. Values remain complex128 and retain
the coefficient unit.

No sparsification, thresholding, conjugation, or unit conversion occurs.

## Numerical validation behavior

None. This record retains a dense coefficient but does not evaluate or compare it.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
