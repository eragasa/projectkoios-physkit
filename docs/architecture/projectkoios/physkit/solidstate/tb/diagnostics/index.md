# Module `projectkoios.physkit.solidstate.tb.diagnostics`

## Current responsibility

The module reports raw scalar or block hopping facts without choosing tolerance or
acceptance policy.

Scalar diagnostics compare coefficients required by Hermiticity or one
caller-supplied integral lattice operation. Block Hermiticity compares
$T[-R]$ with $T[R]^\dagger$ once per available unordered displacement pair and
records a unit-aware Frobenius defect; missing partners remain explicit. Pairing
uses ordinary integer-lattice negation, not equality modulo a finite mesh. An
even-mesh Nyquist representative therefore remains unpaired unless its ordinary
positive partner is present. Block power analysis records the raw squared Frobenius
norm of each hopping matrix with
the squared energy unit and can sum power outside any explicit retained inventory.

## Numerical boundary

Unit-free dense matrix reduction delegates directly to
`projectkoios.numerics.linalg.numpy.frobenius.NumpyFrobeniusNorm` and
`NumpyFrobeniusPower` at the exact Optimization dependency revision. PhysKit
passes immutable coefficient arrays without tuple conversion and restores
hopping-model, basis, displacement, and energy or squared-energy unit semantics.
Nonfinite block subtraction is rejected at its owning arithmetic boundary before
norm evaluation.

## Exclusions

The module does not infer operations, test group closure, add missing partners,
project coefficients, select a displacement inventory, choose a range or cutoff,
normalize power, calculate campaign-specific locality measures, or assign
acceptance.

## Navigation

- [`ScalarHoppingCoefficientComparison`](ScalarHoppingCoefficientComparison/index.md)
- [`ScalarHoppingHermiticityResult`](ScalarHoppingHermiticityResult/index.md)
- [`ScalarHoppingHermiticityAnalyzer`](ScalarHoppingHermiticityAnalyzer/index.md)
- [`ScalarHoppingSymmetryResult`](ScalarHoppingSymmetryResult/index.md)
- [`ScalarHoppingSymmetryAnalyzer`](ScalarHoppingSymmetryAnalyzer/index.md)
- [`BlockHoppingCoefficientComparison`](BlockHoppingCoefficientComparison/index.md)
- [`BlockHoppingHermiticityResult`](BlockHoppingHermiticityResult/index.md)
- [`BlockHoppingHermiticityAnalyzer`](BlockHoppingHermiticityAnalyzer/index.md)
- [`BlockHoppingPower`](BlockHoppingPower/index.md)
- [`BlockHoppingPowerResult`](BlockHoppingPowerResult/index.md)
- [`BlockHoppingPowerAnalyzer`](BlockHoppingPowerAnalyzer/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/`

## Evidence

Tests cover available, exact, defective, and missing conjugate pairs; canonical
pairing; unit-aware block powers; arbitrary explicit omitted inventories; supplied
scalar lattice operations; and floating-point failure behavior. They establish
represented software and bounded numerical behavior only.
