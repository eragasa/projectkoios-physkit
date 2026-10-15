# Class `BlockHoppingBlochEvaluation`

## Current contract

`BlockHoppingBlochEvaluation` correlates the exact `BlockHoppingModel`, original
unitless fractional wave-vector matrix, and one `ComplexMatrixQuantity` per
coordinate row. Every matrix has shape `(model.basis_size, model.basis_size)` and
the model energy unit. Empty coordinate matrices correlate with an empty matrix
inventory.

The record does not claim that returned matrices are Hermitian or diagonalized.

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/bloch.py::BlockHoppingBlochEvaluation`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__BlockHoppingBlochEvaluation__init.py`
