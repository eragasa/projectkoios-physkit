# Class `BlockHoppingBlochEvaluator`

This page defines unit-aware Bloch matrices for one block model or a compatible model
family.

## Current contract

The evaluator maps ordered matrix-entry coefficient series to correlated Bloch
matrices without adding spectral or truncation policy.

- `evaluate(...)` returns one `BlockHoppingBlochEvaluation`.
- `evaluate_family(...)` returns one evaluation per model in supplied order.
- Every family result retains the exact shared coordinate object.
- Models share lattice dimension, basis size and identity, energy unit, and energy
  reference.
- Displacement inventories may differ.
- Empty coordinates return every requested model with an empty matrix tuple.

## Accelerated default

The maintained default delegates to one shared-grid Numerics operation.

- The displacement union is sorted lexicographically.
- Coefficient channels are ordered by model, basis row, and basis column.
- Missing model terms become exact complex zeros on union modes.
- Owner output shape `P × (F n²)` reshapes deterministically to `F × P × n × n`.
- One NumPy phase matrix and one dense matrix multiplication serve the full family.
- Dense mode reduction is tolerance-equivalent rather than bitwise identical to the
  ordered scalar backend.

## Compatibility path

Existing nominal scalar evaluator injection remains supported.

- Supplying `finite_series_evaluator` without a batched evaluator selects entrywise
  ordered evaluation.
- Supplying `batched_finite_series_evaluator` selects the shared-grid owner.
- The scalar and batched properties expose the retained nominal collaborators.
- Both paths preserve the public physical convention, shapes, units, and ordering.

## Exclusions

The adapter performs only Bloch assembly.

- No diagonalization or Hermiticity decision is performed.
- No retained inventory is selected.
- No distance, cutoff, error, normalization, or acceptance result is produced.

## Mapping

The implementation and focused tests are maintained together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bloch.py::BlockHoppingBlochEvaluator`
- Constructor tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__BlockHoppingBlochEvaluator__init.py`
- Evaluation tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__BlockHoppingBlochEvaluator__evaluate.py`
