# Class `BlockHoppingBandDiagnostics`

This page defines aggregate diagnostics for a nonempty reciprocal-point batch.

## Current contract

`BlockHoppingBandDiagnostics` contains componentwise maxima for a nonempty
Bloch-matrix batch.

- `maximum_input_hermiticity_absolute_defect` is the largest point defect.
- `maximum_representative_adjustment_absolute_defect` is the largest point
  adjustment.
- `maximum_input_matrix_absolute_entry` is the largest point input scale.
- `maximum_applied_hermiticity_cutoff` is the largest point cutoff.
- Every field is nonnegative and carries the model energy unit.
- The maximum defect does not exceed the maximum cutoff in a successful result.
- Empty band evaluations use `diagnostics=None` because empty maxima are undefined.

## Mapping

The implementation and invariant tests are maintained together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bands.py::BlockHoppingBandDiagnostics`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingBandDiagnostics__init.py`
