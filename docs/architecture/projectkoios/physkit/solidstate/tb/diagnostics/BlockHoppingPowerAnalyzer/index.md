# Class `BlockHoppingPowerAnalyzer`

## Current contract

`execute(model)` calculates one raw unit-aware Frobenius power for every ordered
block hopping term and returns their correlated `BlockHoppingPowerResult`.
Each immutable dense coefficient array is passed directly to
`NumpyFrobeniusPower(matrix).result()` without tuple conversion. Nonfinite power
arithmetic raises `FloatingPointError`.

The analyzer does not normalize values, select terms, choose a threshold, or
interpret locality.

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/diagnostics.py::BlockHoppingPowerAnalyzer`
- Tests: `tests/projectkoios/physkit/solidstate/tb/diagnostics/test__BlockHoppingPowerAnalyzer__execute.py`
