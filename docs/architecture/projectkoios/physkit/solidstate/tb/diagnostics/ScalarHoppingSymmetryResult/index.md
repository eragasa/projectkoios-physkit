# Class `projectkoios.physkit.solidstate.tb.diagnostics.ScalarHoppingSymmetryResult`

## Current contract

The immutable result retains the exact scalar hopping model, one exact supplied
`IntegralLatticeOperation`, and one comparison per model term in canonical
source order.

`missing_targets` returns absent transformed displacements.
`maximum_absolute_coefficient_difference` returns the largest available
comparison or `None` when every transformed target is missing. The result does
not claim that the operation is a physical or metric symmetry.

## Navigation

- [Parent module](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::ScalarHoppingSymmetryResult`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__ScalarHoppingSymmetryResult__init.py` and `test__ScalarHoppingSymmetryAnalyzer__execute.py`
