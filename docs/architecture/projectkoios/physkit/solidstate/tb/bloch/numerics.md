# Block-hopping Bloch numerics

For fractional reciprocal coordinates $\mathbf q$ and integer displacements
$\mathbf R$, both routes evaluate

$$
H(\mathbf q)=\sum_{\mathbf R}T(\mathbf R)
\exp(+2\pi i\,\mathbf q\cdot\mathbf R).
$$

The exponent is positive and no normalization factor is applied. Dense evaluation
uses either exact ordered entrywise finite-series evaluation or the configured
shared-grid batched backend. Sparse evaluation accumulates phase-weighted CSR blocks,
canonicalizes the CSR result, and rejects nonfinite arithmetic. Results retain the
input coordinate record and model energy unit.

No route diagonalizes the matrix, checks Hermiticity, or applies a tolerance.

## Numerical validation behavior

Numerical tests compare scalar ordered and dense batched evaluation with independently
constructed Fourier sums. Dense/sparse parity tests use equivalent models at shared
points, forbid sparse densification, and verify fixed-unitary basis covariance.
Explicit binary64 tolerances are stated by each comparison; these controls do not
establish material-model adequacy.

Specific pytest evidence:

- `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseBlockHoppingBlochEvaluator__evaluate.py::test__DenseBlockHoppingBlochEvaluator__phase__matches_direct_sum`
- `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseSparseBlockHoppingBlochParity.py::TestDenseSparseBlockHoppingBlochParity::test__SparseBlockHoppingBlochEvaluator__storage__does_not_densify`
- `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseSparseBlockHoppingBlochParity.py::TestDenseSparseBlockHoppingBlochParity::test__DenseSparseBlochEvaluators__representation_parity__agree`
- `tests/projectkoios/physkit/solidstate/tb/bloch/test__DenseSparseBlockHoppingBlochParity.py::TestDenseSparseBlockHoppingBlochParity::test__DenseBlockHoppingBlochEvaluator__covariance__commutes_with_unitary`

## Navigation

- [Module index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
