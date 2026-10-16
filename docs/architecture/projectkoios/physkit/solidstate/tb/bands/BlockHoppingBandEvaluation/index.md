# Class `BlockHoppingBandEvaluation`

This page defines the immutable correlated result returned by band evaluation.

## Current contract

`BlockHoppingBandEvaluation` retains one authoritative Bloch evaluation and every
correlated output of the complete-spectrum calculation.

- `bloch_evaluation` owns the exact model, coordinates, and Bloch matrices.
- `model` and `fractional_wave_vectors` are derived properties rather than copied
  fields.
- `band_energies` has shape `(point_count, basis_size)` and the model energy unit.
- Every energy row is complete and ascending.
- `point_diagnostics[p]` corresponds to coordinate row and Bloch matrix `p`.
- Nonempty `diagnostics` equals the componentwise maxima of point diagnostics.
- Empty evaluations contain shape `(0, basis_size)`, an empty point tuple, and
  `diagnostics=None`.
- Backend strings declare the selected solver and version.

## Interpretation boundary

The record reports a completed reusable numerical operation without adding
consumer-specific scientific meaning.

- No band gap or reference-band selection is present.
- No fitting, error, convergence, or acceptance result is present.
- Structural correlation is validated without a second eigensolver call.
- Numerical derivation remains owned by `BlockHoppingBandEvaluator`.

## Mapping

The implementation and invariant tests are maintained together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bands.py::BlockHoppingBandEvaluation`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingBandEvaluation__init.py`
