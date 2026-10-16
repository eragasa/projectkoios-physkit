# Class `BlockHoppingBandEvaluator`

This page defines unit-aware complete-band evaluation for one model or a compatible
ordered family.

## Current contract

The evaluator adapts block Bloch matrices to complete ascending energy bands under a
caller-supplied numerical tolerance.

- `evaluate_bloch(...)` retains one exact supplied Bloch evaluation.
- `evaluate(...)` evaluates one model's Bloch matrices once.
- `evaluate_bloch_family(...)` concatenates compatible existing Bloch results.
- `evaluate_family(...)` requests one shared Bloch family before diagonalization.
- Family methods preserve model order and exact shared-coordinate identity.
- A strictly excessive defect rejects the whole matrix batch before diagonalization.
- Accepted matrices use $(A+A^\dagger)/2$ and satisfy exact Hermiticity.
- Every accepted nonempty family uses one complete-spectrum batch call.
- Empty shared coordinates return all model results without invoking linear algebra.

## Compatibility requirements

One family represents models in a common physical and numerical space.

- Lattice dimensions agree.
- Ordered-basis sizes and identifiers agree.
- Energy units and energy-reference identifiers agree.
- Every Bloch result retains the exact same coordinate record.
- One tolerance object applies independently to every matrix.

## Failure semantics

Failure classes distinguish malformed inputs, rejected numerical inputs, arithmetic
failures, and backend nonconvergence.

- `TypeError` identifies wrong nominal ownership.
- `ValueError` identifies incompatible family metadata, units, correlation, or a
  defect above its cutoff.
- `FloatingPointError` identifies nonfinite arithmetic or a failed exact-Hermitian
  postcondition.
- `RuntimeError` identifies complete-spectrum solver nonconvergence.
- Exact exception messages remain implementation details.

## Numerical backend

The implementation delegates directly to the pinned canonical Numerics complete-
spectrum contract.

- Evaluation-major concatenation gives each model one contiguous point block.
- `BatchedHermitianSpectrumDefinition.empty` retains basis size for empty points.
- `BatchedHermitianSpectrumTolerance` receives unit-stripped values.
- `NumpyBatchedHermitianSpectrumSolver` is the selected nominal solver.
- The backend identifier is `numpy.linalg.eigvalsh`.
- Per-model energies and diagnostics are partitioned deterministically.
- PhysKit restores physical units and exact Bloch-result correlation.

## Mapping

The implementation and focused construction/evaluation tests are maintained
together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bands.py::BlockHoppingBandEvaluator`
- Tests:
  - `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingBandEvaluator__init.py`
  - `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingBandEvaluator__evaluate.py`
