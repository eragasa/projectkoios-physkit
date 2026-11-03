# `SparseBlockHoppingModel` numerics

The model accepts only exact sparse terms. The abstract owner validates common CSR
block shape, unit, dimension, sorted unique displacements, and identities.
`basis_size` is inherited and derived from the first term.

The model does not densify coefficients, assemble a finite operator, evaluate Bloch
matrices, or choose sparse eigensolver parameters.

## Numerical validation behavior

None. The model retains sparse coefficients but performs no Fourier evaluation,
diagonalization, or numerical comparison.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
