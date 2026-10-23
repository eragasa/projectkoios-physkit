# `RealSymmetricEigendecomposer`

## Responsibility

`RealSymmetricEigendecomposer` strips units from one represented exact
real-symmetric operator, delegates complete or lowest-state execution to the
unit-free `RealSymmetricEigenpairSolver`, and restores the operator unit in one
`RealSymmetricEigendecomposition`.

`execute(operator)` returns the complete spectrum. `execute_lowest(operator,
count)` returns the lowest requested subset and requires `count` to be a
positive built-in `int` no larger than the operator dimension.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/eigenproblems/real_symmetric.py::RealSymmetricEigendecomposer`
- Complete tests: `tests/projectkoios/physkit/linalg/eigenproblems/real_symmetric/test__RealSymmetricEigendecomposer__execute.py`
- Partial tests: `tests/projectkoios/physkit/linalg/eigenproblems/real_symmetric/test__RealSymmetricEigendecomposer__execute_lowest.py`

## Navigation

- [Module](../index.md)
