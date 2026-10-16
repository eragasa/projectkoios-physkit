# Class `BlockHoppingRangeBandEvaluation`

This page defines the correlated output of one source-and-truncation family
calculation.

## Current contract

The result retains full source bands and ordered exact truncation/band pairs.

- `source_band_evaluation` includes the full source Bloch matrices required for later
  model-aware transformations.
- `entries` is nonempty and preserves caller selection order.
- Every truncation retains the exact source model.
- Every band result retains the exact shared coordinate and tolerance objects.
- Retained-model identifiers are unique.
- Retained displacement inventories are unique.
- `source`, `fractional_wave_vectors`, and `tolerance` are derived properties.
- Empty reciprocal coordinates preserve the full result structure with `(0, n)`
  energy arrays and absent aggregate diagnostics.

## Interpretation boundary

The result exposes reusable computations rather than a truncation-quality decision.

- No comparison between source and retained bands is performed.
- No omitted-power normalization is applied.
- No range ranking, stopping rule, or accepted cutoff is selected.
- Consumer-specific error metrics remain outside PhysKit.

## Mapping

The implementation and invariant tests are maintained together.

- Code: `src/python/projectkoios/physkit/solidstate/tb/bands.py::BlockHoppingRangeBandEvaluation`
- Tests: `tests/projectkoios/physkit/solidstate/tb/bands/test__BlockHoppingRangeBandEvaluation__init.py`
