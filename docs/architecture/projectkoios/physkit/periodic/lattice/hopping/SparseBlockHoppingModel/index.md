# Class `SparseBlockHoppingModel`

## Current contract

`SparseBlockHoppingModel` retains a nonempty ordered inventory of
`SparseBlockHoppingTerm` records with one lattice dimension, common square sparse
shape, exact energy unit, energy reference, and ordered-basis identity.
Displacements are unique and lexicographically sorted. `basis_size` is derived from
the common sparse shape.

The model is structural. It neither converts coefficients to dense matrices nor
infers Hermiticity, selects truncations, diagonalizes operators, or declares
scientific acceptance.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py::SparseBlockHoppingModel`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/test__SparseBlockHoppingModel__init.py`
