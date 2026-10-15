# Module `projectkoios.physkit.periodic.lattice.hopping`

## Current responsibility

The module owns immutable translation-invariant scalar and matrix-valued lattice
hopping records plus localized scalar perturbation records. Hopping models retain
model, lattice-dimension, energy-reference, ordered-basis, displacement-order, and
energy-unit identity without imposing Hermiticity or symmetry.

`BlockHoppingTerm` carries one immutable finite nonempty square
`ComplexMatrixQuantity`. `BlockHoppingModel` requires nonempty terms in unique
lexicographic displacement order, one lattice dimension, one square matrix shape,
and one exact energy unit. `basis_size` is derived from that common shape rather
than supplied as redundant state. Displacements are ordinary integer-lattice
translations; the model does not declare a cyclic modulus or finite-DFT alias
classes.

The module does not choose a basis, infer missing partners, transform reciprocal
samples, truncate terms, choose tolerances, or assign physical acceptance.

## Navigation

- [`BlockHoppingTerm`](BlockHoppingTerm/index.md)
- [`BlockHoppingModel`](BlockHoppingModel/index.md)
- [Parent package](../index.md)
- [Source mapping](../../../../../../provenance/finite-periodic-lattice-source-mapping.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/`

## Evidence boundary

Tests establish record invariants, immutable matrix storage, canonical displacement
order, homogeneous shapes and units, and basis-size derivation. They do not validate
a material Hamiltonian or establish that a model is Hermitian.
