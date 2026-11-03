# `SparseBlockHoppingBlochEvaluator` numerics

For each coordinate row, the evaluator starts with an empty complex128 CSR array and
adds every coefficient CSR block multiplied by its complex phase. The result is
canonicalized by summing duplicates, eliminating exact zeros, and sorting indices.
Nonfinite real or imaginary data raises `FloatingPointError`.

The evaluator does not construct a dense matrix, call the dense evaluator, apply a
sparsity threshold, or provide a dense fallback.

## Numerical validation behavior

Sparse evaluation is compared with the dense route for mathematically identical
models over a shared reciprocal grid. The test replaces coefficient `to_dense()` with
a failing sentinel, verifies canonical CSR output, and compares values under an
explicit binary64 tolerance.

Specific pytest evidence:

- `tests/projectkoios/physkit/solidstate/tb/bloch/test__SparseBlockHoppingBlochEvaluator__evaluate.py::test__SparseBlockHoppingBlochEvaluator__gamma__matches_analytic_matrix`
- `tests/projectkoios/physkit/solidstate/tb/bloch/test__SparseBlockHoppingBlochEvaluator__evaluate.py::test__SparseBlockHoppingBlochEvaluator__quarter_turn__matches_matrix`
- `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseSparseBlockHoppingBlochParity.py::TestDenseSparseBlockHoppingBlochParity::test__SparseBlockHoppingBlochEvaluator__storage__does_not_densify`
- `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseSparseBlockHoppingBlochParity.py::TestDenseSparseBlockHoppingBlochParity::test__DenseSparseBlochEvaluators__representation_parity__agree`
- `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseSparseBlockHoppingBlochParity.py::TestDenseSparseBlockHoppingBlochParity::test__SparseBlockHoppingBlochEvaluator__result__is_canonical`

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
