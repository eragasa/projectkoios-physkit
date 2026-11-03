# `DenseBlockHoppingModel` numerics

The model accepts only exact dense terms. The abstract owner validates common shape,
unit, dimension, sorted unique displacements, and identities. `basis_size` is
inherited and derived from the first term.

The model stores no assembled finite operator and performs no Fourier evaluation,
diagonalization, truncation, or coefficient conversion.

## Numerical validation behavior

None. The model retains dense coefficients but performs no Fourier evaluation,
diagonalization, or numerical comparison.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
