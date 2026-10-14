# Class `projectkoios.physkit.solidstate.tb.diagnostics.ScalarHoppingHermiticityAnalyzer`

## Current contract

`execute(model)` compares `t[-R]` with `conjugate(t[R])`. Each available
unordered pair is compared once, including the self-paired zero displacement.
A one-sided hopping produces a comparison whose target is explicit and whose
absolute difference is `None`.

The analyzer reports numerical facts only. It neither adds missing partners nor
chooses a tolerance or Hermiticity status.

## Navigation

- [Parent module](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::ScalarHoppingHermiticityAnalyzer`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__ScalarHoppingHermiticityAnalyzer__execute.py`
