# `DenseBlockHoppingBlochEvaluator` numerics

The ordered route evaluates one scalar finite Fourier series per matrix entry. The
batched route forms a sorted displacement union, inserts exact complex zeros for
missing family terms, and evaluates channels ordered by model, row, and column.
Output reshaping preserves model and coordinate order.

The dense batched reduction is numerically equivalent within floating-point
roundoff but is not required to be bitwise identical to ordered scalar reduction.
Empty coordinate grids produce empty correlated result inventories.

## Numerical validation behavior

Tests compare dense evaluation with an independently assembled finite Fourier sum and
verify that conjugate real-space support remains Hermitian without evaluator-side
symmetrization. A fixed-unitary covariance test verifies row/column and phase
conventions through $H'(\mathbf q)=U^\dagger H(\mathbf q)U$ under stated binary64
tolerances.

Specific pytest evidence:

- `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseBlockHoppingBlochEvaluator__evaluate.py::test__DenseBlockHoppingBlochEvaluator__phase__matches_direct_sum`
- `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseBlockHoppingBlochEvaluator__evaluate.py::test__DenseBlockHoppingBlochEvaluator__support__preserves_hermiticity`
- `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseSparseBlockHoppingBlochParity.py::TestDenseSparseBlockHoppingBlochParity::test__DenseBlockHoppingBlochEvaluator__covariance__commutes_with_unitary`

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
