# Module `projectkoios.physkit.periodic.lattice.hopping`

## Current responsibility

The module owns immutable translation-invariant scalar and matrix-valued lattice
hopping records plus localized scalar perturbation records. Hopping models retain
model, lattice-dimension, energy-reference, ordered-basis, displacement-order, and
energy-unit identity without imposing Hermiticity or symmetry.

`BaseBlockHoppingTerm` and `BaseBlockHoppingModel` are nominal abstract owners for
shared block-hopping identity. `BlockHoppingTerm` and `BlockHoppingModel` are
object-oriented factories that select an exact concrete representation. Dense and
sparse storage remain separate concrete values rather than one coefficient union.

`DenseBlockHoppingTerm` carries one immutable finite nonempty square
`ComplexMatrixQuantity`; `SparseBlockHoppingTerm` carries one canonical immutable
`ComplexSparseMatrixQuantity`. Their model siblings require nonempty terms in unique
lexicographic displacement order, one lattice dimension, one square block shape,
and one exact energy unit. `basis_size` is derived from that common shape rather
than supplied as redundant state. Neither representation silently converts to the
other. Displacements are ordinary integer-lattice translations; the model does not
declare a cyclic modulus or finite-DFT alias classes.

The module does not choose a basis, infer missing partners, transform reciprocal
samples, truncate terms, choose tolerances, or assign physical acceptance.

## Navigation

- [`BaseBlockHoppingTerm`](BaseBlockHoppingTerm/index.md)
- [`BlockHoppingTerm`](BlockHoppingTerm/index.md)
- [`DenseBlockHoppingTerm`](DenseBlockHoppingTerm/index.md)
- [`SparseBlockHoppingTerm`](SparseBlockHoppingTerm/index.md)
- [`BaseBlockHoppingModel`](BaseBlockHoppingModel/index.md)
- [`BlockHoppingModel`](BlockHoppingModel/index.md)
- [`DenseBlockHoppingModel`](DenseBlockHoppingModel/index.md)
- [`SparseBlockHoppingModel`](SparseBlockHoppingModel/index.md)
- [Parent package](../index.md)
- [Source mapping](../../../../../../provenance/finite-periodic-lattice-source-mapping.md)

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Local mapping

- Code: `src/python/projectkoios/physkit/periodic/lattice/hopping.py`
- Tests: `tests/projectkoios/physkit/periodic/lattice/hopping/`

## Evidence boundary

Tests establish record invariants, immutable matrix storage, canonical displacement
order, homogeneous shapes and units, and basis-size derivation. They do not validate
a material Hamiltonian or establish that a model is Hermitian.
