# Class `DenseBlockHoppingBlochEvaluation`

## Current contract

`DenseBlockHoppingBlochEvaluation` correlates the exact `DenseBlockHoppingModel`, original
unitless fractional wave-vector matrix, and one `ComplexMatrixQuantity` per
coordinate row. Every matrix has shape `(model.basis_size, model.basis_size)` and
the model energy unit. Empty coordinate matrices correlate with an empty matrix
inventory.

The record does not claim that returned matrices are Hermitian or diagonalized.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/bloch.py::DenseBlockHoppingBlochEvaluation`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseBlockHoppingBlochEvaluation__init.py`
