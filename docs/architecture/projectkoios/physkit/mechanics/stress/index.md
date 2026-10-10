# Module `projectkoios.physkit.mechanics.stress`

## Responsibility

This module owns the immutable, unit-aware three-dimensional Cartesian stress
quantity used by reusable physics models and outward simulation contracts.

## Public contracts

- `StressTensor` retains one finite 3×3 `MatrixQuantity` whose unit has pressure
  dimensions.
- `StressSignConvention` distinguishes tension-positive and
  compression-positive component conventions.
- Construction preserves every supplied component. It does not symmetrize,
  rotate, convert, or reverse signs implicitly.
- `StressTensor.to_unit()` performs explicit pressure-unit conversion.
- `StressTensor.to_si()` returns the equivalent tensor in pascals.
- `StressTensor.to_sign_convention()` explicitly reverses every component when
  changing convention.
- `StressTensor.hydrostatic_pressure` returns mean normal pressure using the
  conventional positive-under-compression interpretation.

`StressTensor` and its nested `MatrixQuantity` are operationally immutable. The
matrix owns copied, byte-backed, read-only binary64 storage.

## Boundaries

Calculator-native output parsers may retain native component units and sign
conventions. Normalization into another convention must be explicit and must not
discard the native evidence. Serialization and comparison policies belong to
the consuming package rather than this reusable quantity.

## Local mapping

- Code: `src/python/projectkoios/physkit/mechanics/stress.py`
- Tests: `tests/projectkoios/physkit/mechanics/stress/`
- Parent package: [`projectkoios.physkit.mechanics`](../index.md)
