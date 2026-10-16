# Class `BlockHoppingRangeBandEntry`

This page defines one correlated retained-model member of a range-family result.

## Current contract

The entry binds an exact truncation result to the complete band evaluation of its
exact retained model.

- `truncation` retains source identity, selected source terms, and raw omitted power.
- `band_evaluation` retains matrices, energies, diagnostics, units, and backend
  evidence.
- `band_evaluation.model` is the exact `truncation.truncated` object.
- The entry adds no error metric, normalization, convergence result, or acceptance
  interpretation.

## Mapping

The implementation and invariant tests are maintained together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bands.py::BlockHoppingRangeBandEntry`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingRangeBandEntry__init.py`
