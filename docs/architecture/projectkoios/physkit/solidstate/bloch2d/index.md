# `projectkoios.physkit.solidstate.bloch2d`

This package owns two-dimensional Bloch physical models and their explicit
reciprocal-space representations.

## Navigation

- [`cell`](cell/index.md) owns `BlochCell2D`, the unit-aware or explicitly
  unitless direct-cell geometry shared by these models. A displayed model cell
  is not necessarily a primitive Bravais cell.
- [`checkerboard`](checkerboard/index.md) owns the bare-electron plane-wave
  Hamiltonian and band solver for a prescribed square checkerboard potential.

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/bloch2d/`
- Tests: `tests/projectkoios/physkit/solidstate/bloch2d/`
