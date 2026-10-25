# Class `BlockHoppingHermiticityAnalyzer`

## Current contract

`execute(model)` compares $T[-R]$ with $T[R]^\dagger$ in canonical model order.
Available unordered pairs are compared once and missing partners remain explicit.
The analyzer reports energy-unit-bearing Frobenius defects and does not repair,
project, tolerize, or accept the model. Negation is on the ordinary integer lattice;
there is no cyclic modulus or Nyquist self-pairing. Finite-mesh alias Hermiticity
must be checked by a separate cyclic-domain consumer. Each finite dense defect
matrix is passed directly to `NumpyFrobeniusNorm(matrix).result()`.

Nonfinite subtraction is rejected immediately at the subtraction boundary, and
nonfinite norm arithmetic raises `FloatingPointError`; neither produces a
diagnostic value.

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::BlockHoppingHermiticityAnalyzer`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__BlockHoppingHermiticityAnalyzer__execute.py`
