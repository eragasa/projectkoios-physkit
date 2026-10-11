# Module `projectkoios.physkit.solidstate.bloch2d.cell`

## Current responsibility

This module owns `BlochCell2D`, an immutable two-dimensional direct-cell record
for Bloch representations. The two ordered direct vectors are the columns of a
`2×2` `MatrixQuantity`. The class derives cell area, reciprocal basis, and
covariant metric while retaining the selected representation convention.

The word *cell* is intentional. The record does not claim that the displayed
cell is primitive. A model such as the checkerboard potential may have an
additional translation symmetry inside the displayed cell.

## Public contract

`BlochCell2D` supports exactly two conventions:

- `UnitSystem.METAL`: direct-basis magnitudes are in ångströms and derived
  quantities retain the corresponding physical units;
- `UnitSystem.NONDIMENSIONAL`: the direct basis and every derived quantity use
  the first-class `Unitless` record.

`UnitSystem.SI` is outside this contract and is rejected. The direct basis must
have shape `(2, 2)` and pass the scale-invariant numerical-rank criterion
described on the [numeric page](numeric.md). The checkerboard consumer builds a
canonical sparse Hamiltonian by default and selects a sparse eigensolver when a
proper partial spectrum is requested.

## Architecture views

- [Schematics](schematics.md): ownership, validation, derivation, and consumer
  relationships.
- [Numeric](numeric.md): rank criterion, dual basis, metric, area, square-cell
  normalization, and nondimensional kinetic convention.
- [Scientific](scientific.md): interpretation, assumptions, limitations, and
  evidence boundaries.

## Navigation

- [Checkerboard model](../checkerboard/index.md)
- [Parent package](../index.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/bloch2d/cell.py`
- Tests: `tests/projectkoios/physkit/solidstate/bloch2d/cell/`

## Evidence boundary

Mapped tests establish the declared type, unit, shape, numerical-rank, and
geometric identities for selected cases. They do not establish that a supplied
cell is primitive, that it represents a real material, or that a downstream
basis truncation is converged.
