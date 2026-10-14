# Module `projectkoios.physkit.solidstate.tb.diagnostics`

## Current responsibility

The module directly compares scalar hopping coefficients required by
Hermiticity or by one caller-supplied integral lattice operation.

## Public contract

A comparison names its source and target displacements. Its absolute coefficient
difference is a finite nonnegative float when the target exists and `None` when
the target is missing. Hermiticity uses
`|t[-R] - conjugate(t[R])|`; supplied-operation symmetry uses
`|t[M R] - t[R]|`.

Each analyzer returns one immutable result retaining the model, its ordered
comparisons, and the supplied operation where applicable. Results expose missing
targets and an optional maximum available difference. The module does not infer
operations, test metric symmetry or group closure, choose tolerances, project
coefficients, or assign acceptance.

## Navigation

- [`ScalarHoppingCoefficientComparison`](ScalarHoppingCoefficientComparison/index.md)
- [`ScalarHoppingHermiticityResult`](ScalarHoppingHermiticityResult/index.md)
- [`ScalarHoppingHermiticityAnalyzer`](ScalarHoppingHermiticityAnalyzer/index.md)
- [`ScalarHoppingSymmetryResult`](ScalarHoppingSymmetryResult/index.md)
- [`ScalarHoppingSymmetryAnalyzer`](ScalarHoppingSymmetryAnalyzer/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/`

## Evidence

Tests cover available, zero, nonzero, and missing comparisons; canonical
Hermitian pairing; supplied signed-permutation and general unimodular
operations; dimensional validation; and optional maximum semantics. No test
establishes scientific acceptance or physical symmetry.
