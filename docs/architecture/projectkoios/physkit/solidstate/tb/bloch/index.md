# Module `projectkoios.physkit.solidstate.tb.bloch`

## Current responsibility

The module owns the physics-semantic bridge from an ordered
`ScalarHoppingModel` to scalar Bloch-Hamiltonian values. For dimensionless
fractional reciprocal coordinates $\mathbf q$ and integer direct-lattice
displacements $\mathbf R$, it evaluates

$$
h(\mathbf q)=\sum_{\mathbf R} t_{\mathbf R}
\exp\!\left(+2\pi i\,\mathbf q\cdot\mathbf R\right).
$$

The positive exponent, factor of $2\pi$, supplied hopping-term order, and absence
of normalization are fixed conventions. The model owns coefficient order,
energy unit, energy-reference identity, and basis identity. The supplied wave
vectors are rows of a unitless `MatrixQuantity` with width equal to the model
dimension.

The result remains complex because `ScalarHoppingModel` does not impose
Hermiticity. Diagnostics, tolerances, acceptance, and interpretation of a small
imaginary component remain separate responsibilities.

## Numerics boundary

`ScalarHoppingBlochEvaluator` converts the model and wave vectors to a
`FiniteFourierSeriesEvaluationDefinition`, calls a nominal
`FiniteFourierSeriesEvaluator`, and wraps the returned built-in complex values in
the model energy unit. The dependency is bound to exact Numerics owner commit
`d9a82252faa1b41d82cce11149329aa91b84dff0`.

The upstream Numerics operation remains unit-free and does not receive lattice,
hopping, Bloch, basis, energy-reference, or acceptance meaning. PhysKit is not a
pass-through facade: this adapter exists because it restores those physical
semantics around the generic finite-series calculation.

The default adapter constructs the reviewed ordered NumPy backend through the
Numerics factory. A nominal evaluator may be injected for verification or a
separately reviewed implementation. Empty wave-vector matrices of shape
`(0, dimension)` produce an empty complex vector carrying the model energy unit.
Numerics validation and numerical failures propagate without broad translation.

## Public classes

- [`ScalarHoppingBlochEvaluator`](ScalarHoppingBlochEvaluator/index.md)
- [`ScalarHoppingBlochEvaluation`](ScalarHoppingBlochEvaluation/index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/bloch.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bloch/`
