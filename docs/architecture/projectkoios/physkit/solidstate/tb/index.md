# Package `projectkoios.physkit.solidstate.tb`

## Current responsibility

The package owns reusable scalar and matrix-valued tight-binding Bloch evaluation,
raw coefficient diagnostics, unit-aware block powers, and explicit
conjugate-closed block-hopping selection. Inputs retain lattice displacement,
ordered basis, model, energy-unit, and energy-reference identity.

Its APIs follow represented equations directly. Callers supply models,
coordinates, operations, and retained displacement inventories to small numerical
actions; no workflow request, campaign state, inferred symmetry family, cutoff,
tolerance, normalized locality measure, or scientific acceptance policy is owned
here.

Unit-free finite Fourier-series and complex-matrix Frobenius kernels belong to
`projectkoios.numerics`. PhysKit delegates to those pinned owner APIs and supplies
lattice, hopping, basis, and unit semantics around their calculations.

## Navigation

- [Module `bloch`](bloch/index.md)
- [Module `diagnostics`](diagnostics/index.md)
- [Module `truncation`](truncation/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/`
- Tests: `tests/projectkoios/physkit/solidstate/tb/`
