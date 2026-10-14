# Class `ScalarHoppingBlochEvaluator`

## Ownership

`ScalarHoppingBlochEvaluator` is the PhysKit adapter for evaluating one
`ScalarHoppingModel` at supplied fractional reciprocal wave vectors. Its public
`evaluate(model, fractional_wave_vectors)` operation follows the Bloch equation
directly and does not introduce a workflow request record.

The evaluator validates that coordinates are unitless and have one column per
model dimension. It maps coordinate rows, ordered integer displacements, and
ordered complex hopping coefficients to the exact dependency-free Numerics
record. The returned complex values are wrapped in the model energy unit and
correlated with the original model and coordinate quantity.

The constructor accepts a nominal `FiniteFourierSeriesEvaluator`. Omitting it
selects the ordered `NUMPY_DENSE` implementation through the Numerics factory.
The adapter does not catch or reinterpret Numerics programming or floating-point
failures.

## Exclusions

The class does not infer missing Hermitian partners, reduce symmetry orbits,
choose wave-vector paths, diagonalize multi-orbital matrices, select a tolerance,
or decide whether a model or result is physically acceptable.

## Mapping

- Implementation: `src/python/projectkoios/physkit/solidstate/tb/bloch.py`
- Constructor tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__ScalarHoppingBlochEvaluator__init.py`
- Evaluation tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__ScalarHoppingBlochEvaluator__evaluate.py`
