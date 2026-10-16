# Class `BlockHoppingBandEvaluator`

This page defines the unit-aware action from Bloch matrices to complete bands.

## Current contract

`BlockHoppingBandEvaluator` adapts block Bloch matrices to complete unit-aware
energy bands under caller-supplied numerical tolerances.

- `evaluate_bloch(bloch_evaluation, tolerance)` is the primitive action.
- The primitive retains and reuses the exact supplied Bloch evaluation.
- `evaluate(model, fractional_wave_vectors, tolerance)` evaluates Bloch matrices
  once and delegates to the primitive.
- A strictly excessive defect rejects the whole batch before diagonalization.
- Accepted matrices produce the explicit representative $(A+A^\dagger)/2$.
- Every representative satisfies the exact-Hermitian postcondition.
- Every accepted nonempty batch uses one batched complete-spectrum call.
- Empty batches return without invoking linear algebra.

## Failure semantics

Failure classes distinguish malformed inputs, rejected numerical inputs, invalid
floating-point results, and backend nonconvergence.

- `TypeError` identifies wrong nominal input ownership.
- `ValueError` identifies incompatible units or a defect above its applied cutoff.
- `FloatingPointError` identifies nonfinite arithmetic/output or failure of the
  exact-Hermitian postcondition.
- `RuntimeError` identifies complete-spectrum solver nonconvergence.
- Exact exception messages remain implementation details.

## Numerical backend

The default implementation delegates directly to the pinned canonical Numerics
complete-spectrum contract.

- `BatchedHermitianSpectrumDefinition` represents nonempty batches.
- `BatchedHermitianSpectrumDefinition.empty` represents empty batches with an
  explicit basis size.
- `BatchedHermitianSpectrumTolerance` receives unit-stripped tolerances.
- `NumpyBatchedHermitianSpectrumSolver` is the selected nominal solver.
- The backend identifier is `numpy.linalg.eigvalsh`.
- The backend version is `numpy.__version__`.
- `UPLO="L"` is explicit for nonempty solves.
- PhysKit restores physical units and Bloch-result correlation.

## Mapping

The implementation and focused construction/evaluation tests are maintained
together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bands.py::BlockHoppingBandEvaluator`
- Tests:
  - `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingBandEvaluator__init.py`
  - `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingBandEvaluator__evaluate.py`
