# `ScalarHoppingBlochEvaluator` numerics

The evaluator maps ordered displacements and complex scalar coefficients to the
unit-free finite-series definition and delegates the ordered reduction to a nominal
Numerics evaluator. The default is the reviewed `NUMPY_DENSE` implementation.

Coordinate rows are interpreted in turns relative to primitive reciprocal axes. The
evaluated complex values are restored to the model energy unit without rounding,
real-part projection, or tolerance filtering. Empty coordinate grids remain valid.

## Numerical validation behavior

Evaluator tests compare delegated results with explicit scalar Bloch sums at selected
coordinates, verify the positive-exponent convention without normalization, and
require nonfinite backend failures to propagate. Any comparison tolerance is declared
by the test rather than embedded in the evaluator.

Specific pytest evidence:

- `tests/projectkoios/physkit/solidstate/tb/bloch/test__ScalarHoppingBlochEvaluator__evaluate.py::test__ScalarHoppingBlochEvaluator__phase__matches_direct_sum`
- `tests/projectkoios/physkit/solidstate/tb/bloch/test__ScalarHoppingBlochEvaluator__evaluate.py::test__ScalarHoppingBlochEvaluator__nonfinite_phase__propagates_failure`

- [Index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientifics.md)
