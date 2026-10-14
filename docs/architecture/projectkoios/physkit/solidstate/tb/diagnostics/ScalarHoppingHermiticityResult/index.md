# Class `projectkoios.physkit.solidstate.tb.diagnostics.ScalarHoppingHermiticityResult`

## Current contract

The immutable result retains the exact scalar hopping model and comparisons in
canonical source order. Available unordered pairs `{R, -R}` occur once; missing
targets remain explicit comparisons.

`missing_targets` returns absent negative displacements.
`maximum_absolute_coefficient_difference` returns the largest available
comparison or `None` when every target is missing. Neither property makes a
Hermiticity decision.

## Navigation

- [Parent module](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::ScalarHoppingHermiticityResult`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__ScalarHoppingHermiticityResult__init.py` and `test__ScalarHoppingHermiticityAnalyzer__execute.py`
