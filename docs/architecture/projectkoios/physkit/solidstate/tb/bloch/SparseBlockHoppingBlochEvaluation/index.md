# Class `SparseBlockHoppingBlochEvaluation`

## Current contract

`SparseBlockHoppingBlochEvaluation` correlates one exact sparse model, one unitless
fractional-coordinate record, and one `ComplexSparseMatrixQuantity` per reciprocal
point. Every matrix has the model basis shape and energy unit. Sparse storage is
part of the nominal result contract, not an optimization hidden behind a dense API.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/bloch.py::SparseBlockHoppingBlochEvaluation`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__SparseBlockHoppingBlochEvaluator__evaluate.py`
