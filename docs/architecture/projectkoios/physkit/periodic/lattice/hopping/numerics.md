# Block-hopping numerics

Block-hopping records own validated representations, not a numerical solver. Dense
coefficients use immutable complex128 matrices. Sparse coefficients use canonical
immutable complex128 CSR storage with sorted indices, summed duplicates, and no
explicit zeros.

Both model representations require a nonempty lexicographically ordered unique
displacement inventory, one common nonempty square coefficient shape, and one exact
unit. `basis_size` is derived from that shape. Factories select exact nominal
representations and perform no dense–sparse conversion.

Hermiticity, Fourier evaluation, diagonalization, truncation, fitting, and tolerance
selection belong to separate actions.

## Navigation

- [Module index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
