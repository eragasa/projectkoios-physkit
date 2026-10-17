# Class `projectkoios.physkit.periodic.lattice.bain.path.DiamondDecoratedBainPath3D`

## Current contract

The stateless action carries a species-neutral two-site decoration along the Bain path.

- `evaluate(...)` delegates interval and Bravais geometry to `UnitVolumeBainPath3D`.
- Every result is a new `DiamondDecoratedBainPathPoint3D`.
- The fractional decoration remains fixed in Bain coordinates.
- Cartesian offsets deform affinely with the current direct basis.
- No atomic species or physical lattice parameter is selected.

## Physical scope

The action represents decorated crystallographic geometry only.

- At the FCC endpoint, the coordinate bridge recovers the conventional nonzero diamond
  offset `(1/4, 1/4, 1/4)`.
- The decoration remains mathematically defined through the complete Bain interval.
- The endpoint at parameter one is not assigned a material phase identity.
- Electronic Hamiltonians, hopping response, energetics, and transformation kinetics
  remain outside this action.

## Local mapping

Source and tests remain adjacent to the module owner.

- Code: `src/python/projectkoios/physkit/periodic/lattice/bain/path.py::DiamondDecoratedBainPath3D`
- Tests: `tests/projectkoios/physkit/periodic/lattice/bain/path/test__DiamondDecoratedBainPath3D__evaluate.py`
- [Parent module](../index.md)
