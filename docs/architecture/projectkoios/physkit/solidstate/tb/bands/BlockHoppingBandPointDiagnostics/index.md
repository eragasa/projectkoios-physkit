# Class `BlockHoppingBandPointDiagnostics`

This page defines the raw diagnostics correlated with one reciprocal point.

## Current contract

`BlockHoppingBandPointDiagnostics` correlates one reciprocal-point index with raw
maximum-absolute numerical diagnostics.

- `point_index` is the canonical zero-based coordinate-row index.
- `input_hermiticity_maximum_absolute_defect` reports
  $\max_{ij}|A_{ij}-A_{ji}^*|$.
- `representative_adjustment_maximum_absolute_defect` reports
  $\max_{ij}|[(A+A^\dagger)/2]_{ij}-A_{ij}|$.
- `input_matrix_maximum_absolute_entry` supplies the relative-tolerance scale.
- `applied_hermiticity_cutoff` records the pointwise numerical cutoff.
- All four quantities are nonnegative and carry one common energy unit.
- Successful results satisfy `defect <= cutoff`.

## Mapping

The implementation and invariant tests are maintained together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bands.py::BlockHoppingBandPointDiagnostics`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingBandPointDiagnostics__init.py`
