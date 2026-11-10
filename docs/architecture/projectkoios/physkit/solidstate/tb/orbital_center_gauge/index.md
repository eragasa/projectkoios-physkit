# Module `projectkoios.physkit.solidstate.tb.orbital_center_gauge`

This module owns the generic diagonal transformation from an ordinary-cell Bloch
representation to a representation carrying explicit fractional orbital-center phases.

## Current contract

- `CellToFullCenterGaugeTransformationDefinition` retains the exact ordered
  `(basis_size, 3)` fractional-center matrix and one three-component fractional wave
  vector.
- Both inputs are immutable nondimensional PhysKit quantities; no raw array boundary is
  exposed.
- `CellToFullCenterGaugeTransformationConstructor` uses the positive phase convention

  $$
  D_{aa}(q)=\exp(+2\pi i q\cdot\tau_a).
  $$

- `CellToFullCenterGaugeTransformation` correlates the exact definition with one
  immutable nondimensional complex diagonal matrix.
- Public records and action inputs require named arguments.

For a cell-gauge matrix, the represented full-center matrix is

$$
H_{\mathrm{full}}(q)=D(q)^\dagger H_{\mathrm{cell}}(q)D(q).
$$

## Class architecture

- [`CellToFullCenterGaugeTransformationDefinition`](CellToFullCenterGaugeTransformationDefinition/index.md)
- [`CellToFullCenterGaugeTransformation`](CellToFullCenterGaugeTransformation/index.md)
- [`CellToFullCenterGaugeTransformationConstructor`](CellToFullCenterGaugeTransformationConstructor/index.md)

## Exclusions

The module does not choose orbital centers, infer basis compatibility, assemble hopping
coefficients, evaluate spectra, define reciprocal sewing policy, or judge scientific
admissibility.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/orbital_center_gauge.py`
- Tests: `tests/projectkoios/physkit/solidstate/tb/orbital_center_gauge/`
- [Parent tight-binding package](../index.md)
