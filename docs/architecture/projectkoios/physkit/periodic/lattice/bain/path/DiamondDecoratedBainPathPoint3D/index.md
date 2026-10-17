# Class `projectkoios.physkit.periodic.lattice.bain.path.DiamondDecoratedBainPathPoint3D`

## Current contract

The immutable point correlates Bravais Bain geometry with three basis representations.

- `bain_point` is an owned `UnitVolumeBainPathPoint3D` copy.
- `fractional_basis_bain` stores the two fixed Bain-coordinate offsets.
- `fractional_basis_fcc` stores their exact image in the declared reference FCC basis.
- `cartesian_basis` stores the current dimensionless Cartesian offsets.
- Every basis is an immutable, explicitly unitless `MatrixQuantity`.

## Correlation

All representations derive from one Bravais path point.

- Direct lattice, metric, and primitive volume delegate to the owned Bain point.
- Cartesian offsets equal the Bain fractional rows multiplied by `A.T`.
- FCC fractional offsets derive from `FccBainCoordinateBridge3D`.
- No independent caller-supplied basis can create inconsistent point state.

## Ownership boundary

The point stops before chemical or electronic modeling.

- Atomic symbols, occupations, pseudopotentials, and physical scale are absent.
- Consumers may compose the fractional geometry with separately owned structure data.
- Synthetic response models, selected bonds, thresholds, and interpretation remain
  consumer responsibilities.

## Local mapping

Source and tests remain adjacent to the module owner.

- Code: `src/python/projectkoios/physkit/periodic/lattice/bain/path.py::DiamondDecoratedBainPathPoint3D`
- Tests: `tests/projectkoios/physkit/periodic/lattice/bain/path/test__DiamondDecoratedBainPath3D__evaluate.py`
- [Parent module](../index.md)
