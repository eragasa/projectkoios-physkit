# Package `projectkoios.physkit.solidstate.tb`

## Current responsibility

The package owns reusable scalar tight-binding Bloch evaluation and diagnostics
built from explicit hopping models, fractional reciprocal wave vectors, and
caller-supplied lattice operations.

Its APIs follow the represented equations directly. Callers pass models,
coordinates, and operations to small numerical actions; no workflow request,
campaign state, symmetry inference, or scientific acceptance policy is owned
here. The unit-free finite Fourier-series kernel belongs to
`projectkoios.numerics`. PhysKit supplies the lattice meaning, hopping-model
identity, and energy unit around that generic calculation.

## Navigation

- [Module `bloch`](bloch/index.md)
- [Module `diagnostics`](diagnostics/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/`
- Tests: `tests/projectkoios/physkit/solidstate/tb/`
