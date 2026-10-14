# Class `projectkoios.physkit.solidstate.tb.diagnostics.ScalarHoppingCoefficientComparison`

## Current contract

The immutable record names a source displacement and its negative or
operation-transformed target. `absolute_coefficient_difference` is a finite
nonnegative float when the target coefficient exists and `None` when it is
missing. Consequently, an exact zero and a missing target remain distinct.

`target_is_present` exposes that distinction without applying a tolerance.

## Navigation

- [Parent module](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::ScalarHoppingCoefficientComparison`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__ScalarHoppingCoefficientComparison__init.py`
