# Class `BlockHoppingTruncationResult`

## Current contract

The result correlates the exact source model, canonical retained source
displacements, distinctly identified truncated model, and raw omitted Frobenius
power. The truncated identifier must differ from the source identifier because the
retained model is a different mathematical object. Dimension, ordered-basis
identity, energy unit, and energy reference are preserved. Truncated terms are the
same immutable term objects owned by the source. `omitted_displacements` exposes
the complementary source inventory.

Construction revalidates conjugate closure and all correlations. No tolerance or
acceptance disposition is stored.

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/truncation.py::BlockHoppingTruncationResult`
- Tests: `tests/projectkoios/physkit/solidstate/tb/truncation/test__BlockHoppingTruncationResult__init.py`
