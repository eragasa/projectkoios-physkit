# Module `projectkoios.physkit.solidstate.tb.bloch`

This page defines the physics-semantic bridge from hopping models to Bloch values.

## Current responsibility

The module evaluates scalar and matrix-valued Bloch Hamiltonians for unitless
fractional reciprocal coordinates $\mathbf q$ and integer displacements $\mathbf R$.

- Scalar models evaluate
  $h(\mathbf q)=\sum_{\mathbf R}t_{\mathbf R}\exp(+2\pi i\,\mathbf q\cdot\mathbf R)$.
- Block models apply the same finite series entrywise to $T_{\mathbf R}$.
- The positive exponent and factor of $2\pi$ are fixed conventions.
- Hopping-term order is preserved and no normalization is applied.
- Models retain coefficient order, energy unit, energy reference, and basis identity.
- Coordinates are rows of a unitless `MatrixQuantity` whose width is the lattice
  dimension.
- Results remain complex because hopping models do not impose Hermiticity.
- Diagnostics, tolerances, diagonalization, acceptance, and interpretation remain
  separate responsibilities.

## Numerics boundary

The adapters surround the pinned unit-free finite-series owner with lattice,
hopping, basis, and unit semantics.

- The exact owner revision is commit
  `fff068dfe38309db05bf6de075bac53bb6597813`, tree
  `4ec1c3d6334d884e4443a779f40321eef7f76981`.
- The scalar adapter constructs one `FiniteFourierSeriesEvaluationDefinition`.
- The block adapter constructs one definition per matrix entry with shared
  coordinates and displacements.
- Each adapter calls a nominal `FiniteFourierSeriesEvaluator`.
- Returned built-in complex values regain the model energy unit.
- The Numerics owner receives no lattice, hopping, Bloch, basis, energy-reference,
  or acceptance meaning.
- Default adapters use the reviewed ordered NumPy backend through the Numerics
  factory.
- Nominal evaluator injection supports verification or another reviewed backend.
- Empty coordinate shape `(0, dimension)` produces an empty scalar vector or block
  matrix inventory.
- Numerics validation and numerical failures propagate without broad translation.

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
