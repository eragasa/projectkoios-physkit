# Class `projectkoios.physkit.solidstate.tb.diagnostics.ScalarHoppingSymmetryAnalyzer`

## Current contract

`execute(model, operation)` compares every source coefficient `t[R]` with the
coefficient at the explicitly transformed displacement `M R`. The result uses
`|t[M R] - t[R]|` when that target exists and records `None` otherwise.

The caller supplies exactly one `IntegralLatticeOperation` per call. The
analyzer does not build an operation family, infer symmetry, test group closure
or metric compatibility, enumerate orbits, project coefficients, or choose a
tolerance.

## Navigation

- [Parent module](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::ScalarHoppingSymmetryAnalyzer`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__ScalarHoppingSymmetryAnalyzer__execute.py`
