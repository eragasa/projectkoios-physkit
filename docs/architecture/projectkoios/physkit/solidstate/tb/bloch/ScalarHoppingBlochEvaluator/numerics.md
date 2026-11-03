# `ScalarHoppingBlochEvaluator` numerics

The evaluator maps ordered displacements and complex scalar coefficients to the
unit-free finite-series definition and delegates the ordered reduction to a nominal
Numerics evaluator. The default is the reviewed `NUMPY_DENSE` implementation.

Coordinate rows are interpreted in turns relative to primitive reciprocal axes. The
evaluated complex values are restored to the model energy unit without rounding,
real-part projection, or tolerance filtering. Empty coordinate grids remain valid.

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
