# Class `BlockHoppingBlochEvaluator`

## Ownership

`BlockHoppingBlochEvaluator` evaluates one `BlockHoppingModel` at supplied
fractional reciprocal wave vectors. It maps every ordered matrix-entry
coefficient series to the generic Numerics finite-series owner and reconstructs
one immutable energy-unit-bearing Bloch matrix per coordinate row.

The constructor accepts a nominal `FiniteFourierSeriesEvaluator`; omission selects
the ordered `NUMPY_DENSE` backend through the Numerics factory. Coordinates must
be unitless and have one column per model dimension. Matrix-entry evaluation is
row-major and shares exact coordinate and displacement order.

The adapter does not diagonalize matrices, infer Hermiticity, select a basis or
wave-vector path, truncate hoppings, choose tolerances, or assign acceptance.

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/bloch.py::BlockHoppingBlochEvaluator`
- Constructor tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__BlockHoppingBlochEvaluator__init.py`
- Evaluation tests: `tests/projectkoios/physkit/solidstate/tb/bloch/test__BlockHoppingBlochEvaluator__evaluate.py`
