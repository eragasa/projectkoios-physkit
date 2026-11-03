# `ScalarHoppingBlochEvaluation` numerics

The result requires a unitless coordinate matrix with width equal to the model
lattice dimension and a complex vector with exactly one value per coordinate row.
The vector unit must equal the model energy unit. Empty coordinate and value
inventories are correlated explicitly.

Construction validates the correlation but does not recompute the finite Fourier
series or alter complex values.

## Numerical validation behavior

None. This result record retains evaluator output but does not independently recompute
or numerically validate the scalar Bloch sum.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
