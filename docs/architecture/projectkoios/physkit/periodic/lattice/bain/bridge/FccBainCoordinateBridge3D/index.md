# Class `projectkoios.physkit.periodic.lattice.bain.bridge.FccBainCoordinateBridge3D`

## Current contract

The immutable record owns one exact integral FCC--Bain coordinate bridge.

- `bain_to_fcc_matrix` implements `n_F = U n_B`.
- `fcc_to_bain_matrix` is the exact integer inverse.
- Both matrices are immutable signed-64-bit arrays.
- The determinant is exactly one.
- Mapping methods accept individual vectors or batches with final dimension three.

## Crystallographic convention

The bridge fixes one declared pair of primitive-coordinate conventions.

- The Bain-to-FCC matrix has rows `(1, 1, 0)`, `(0, 0, 1)`, and `(0, -1, 0)`.
- The nonzero decorated Bain offset `(1/2, -1/4, 1/4)` maps to the FCC diamond offset
  `(1/4, 1/4, 1/4)`.
- Integer maps round-trip exactly within represented signed-64-bit arithmetic.
- Unimodularity preserves a bijection modulo every finite periodic size.
- Fractional coordinates are not silently reduced modulo one.

## Evidence limits

The map establishes coordinate identity rather than material behavior.

- Tests bind exact matrices, inverse composition, the diamond offset, and quotient
  bijections for multiple finite sizes.
- The bridge does not prove that two independently supplied lattices use these bases.
- It does not attach a material, deformation mechanism, or scientific acceptance.

## Local mapping

Source and focused verification remain adjacent.

- Code: `src/python/projectkoios/physkit/periodic/lattice/bain/bridge.py::FccBainCoordinateBridge3D`
- Tests: `tests/projectkoios/physkit/periodic/lattice/bain/bridge/test__FccBainCoordinateBridge3D__contract.py`
- [Parent module](../index.md)
