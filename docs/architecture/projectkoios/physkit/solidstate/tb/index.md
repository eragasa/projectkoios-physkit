# Package `projectkoios.physkit.solidstate.tb`

This page defines the reusable tight-binding package boundary.

## Current responsibility

The package owns reusable tight-binding representations, numerical adapters, and
raw diagnostics with explicit physics semantics.

- Scalar and matrix-valued hopping models retain lattice-displacement identity.
- Bloch evaluation retains fractional-coordinate and energy-unit conventions.
- Complete block-band evaluation retains ordered-basis and backend evidence.
- Hermiticity and power diagnostics report raw numerical facts.
- Conjugate-closed truncation preserves explicit model and term identity.
- Workflow state, campaign policy, scientific tolerances, cutoffs, normalized
  locality metrics, and acceptance decisions remain outside the package.
- Unit-free finite-series, Frobenius, and batched Hermitian-spectrum kernels belong
  to `projectkoios.numerics`.
- Published pinned Numerics APIs own all integrated generic kernels.
- Complete block-band evaluation delegates directly to the canonical batched
  Hermitian-spectrum definition, tolerance, diagnostics, result, and NumPy solver.

## Navigation

The module pages describe each cohesive tight-binding capability.

- [Module `bloch`](bloch/index.md)
- [Module `bands`](bands/index.md)
- [Module `diagnostics`](diagnostics/index.md)
- [Module `truncation`](truncation/index.md)
- [Parent package](../index.md)

## Local mapping

The source and tests share one package boundary.

- Code: `src/python/projectkoios/physkit/solidstate/tb/`
- Tests: `tests/projectkoios/physkit/solidstate/tb/`
