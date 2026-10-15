# Class `BlockHoppingModel`

## Current contract

`BlockHoppingModel` retains one nonempty ordered inventory of `BlockHoppingTerm`
records together with model identity, lattice dimension, exact energy unit,
energy-reference identity, and ordered-basis identity. Displacements are unique and
lexicographically sorted. Every term has the model dimension, common square shape,
and exact model unit.

`basis_size` is the common coefficient row and column count. It is derived rather
than separately supplied.

The model is structural: it does not require conjugate partners, test coefficient
Hermiticity, select terms, diagonalize Bloch matrices, or encode campaign policy.
Its displacements belong to the ordinary integer lattice. Centered finite-DFT
representatives, cyclic equivalence, and even-mesh Nyquist self-pairing require a
separate consumer-owned cyclic-domain contract.

## Mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py::BlockHoppingModel`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/test__BlockHoppingModel__init.py`
