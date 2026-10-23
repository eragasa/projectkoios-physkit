# `VerifiedMatrixUtils`

## Responsibility

`VerifiedMatrixUtils` provides one uniform explicit dense-materialization
boundary for the verified matrix-quantity union. Dense input is returned by
identity; sparse input is materialized as a new immutable dense quantity with
the exact retained unit.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/matrices.py::VerifiedMatrixUtils`
- Tests: `tests/projectkoios/physkit/linalg/matrices/test__VerifiedMatrixUtils__to_dense_matrix.py`

## Navigation

- [Module](../index.md)
