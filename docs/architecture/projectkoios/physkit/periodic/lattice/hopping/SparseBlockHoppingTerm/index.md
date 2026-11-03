# Class `SparseBlockHoppingTerm`

## Current contract

`SparseBlockHoppingTerm` retains one ordinary-lattice displacement and one exact
`ComplexSparseMatrixQuantity`. The coefficient must have a nonempty square shape.
Its canonical immutable CSR arrays and energy unit remain owned by the PhysKit
quantity; the term does not densify, infer partners, or apply a tolerance.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py::SparseBlockHoppingTerm`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/test__SparseBlockHoppingTerm__init.py`
