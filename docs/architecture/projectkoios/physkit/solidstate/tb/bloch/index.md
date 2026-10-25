# Module `projectkoios.physkit.solidstate.tb.bloch`

This page defines physics-semantic bridges from hopping models to Bloch values and
matrices.

## Current responsibility

The module evaluates scalar and matrix-valued Bloch Hamiltonians for unitless
fractional reciprocal coordinates $\mathbf q$ and integer displacements $\mathbf R$.

- Scalar models evaluate
  $h(\mathbf q)=\sum_{\mathbf R}t_{\mathbf R}\exp(+2\pi i\,\mathbf q\cdot\mathbf R)$.
- Block models apply the same finite series to every ordered matrix entry.
- The positive exponent and absence of normalization are fixed conventions.
- Models retain coefficient order, energy unit, energy reference, and basis identity.
- Results remain complex because hopping models do not impose Hermiticity.
- Diagnostics, tolerances, diagonalization, acceptance, and interpretation remain
  separate responsibilities.

## Scalar Numerics boundary

The scalar adapter preserves the exact ordered finite-series contract.

- It constructs one `FiniteFourierSeriesEvaluationDefinition`.
- The default scalar evaluator uses the reviewed ordered NumPy backend.
- Point, term, and coordinate-component reduction order remain observable.
- Empty coordinate shape `(0, dimension)` produces an empty complex vector.

## Block and family Numerics boundary

The default block adapter uses the pinned shared-grid batched finite-Fourier owner.

- The exact owner revision is commit
  `31ce3de66a46302e21b87189b1af4ce78adc0f7d`, tree
  `44b33434f4a2226045ef9db35704567df4b78b6f`.
- `evaluate(...)` maps one model to `basis_size²` coefficient channels.
- `evaluate_family(...)` maps models in supplied order to one combined channel set.
- Compatible families share dimension, ordered basis, energy unit, and energy
  reference.
- The sorted displacement union is the shared integer-mode inventory.
- Exact complex zeros represent modes absent from one family member.
- `NumpyBatchedFiniteFourierSeriesEvaluator` forms one point-mode phase matrix and
  performs one dense channel reduction.
- Dense reduction may differ from the scalar ordered path by binary64 rounding.
- Supplying only a scalar evaluator retains the previous entrywise compatibility
  path.
- PhysKit restores each model's exact identity, coordinate object, matrix shape, and
  energy unit.

## Public classes

The public records and evaluators distinguish scalar and block ownership.

- [`ScalarHoppingBlochEvaluator`](ScalarHoppingBlochEvaluator/index.md)
- [`ScalarHoppingBlochEvaluation`](ScalarHoppingBlochEvaluation/index.md)
- [`BlockHoppingBlochEvaluator`](BlockHoppingBlochEvaluator/index.md)
- [`BlockHoppingBlochEvaluation`](BlockHoppingBlochEvaluation/index.md)

## Local mapping

The implementation and mapped tests share one module boundary.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bloch.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bloch/`
