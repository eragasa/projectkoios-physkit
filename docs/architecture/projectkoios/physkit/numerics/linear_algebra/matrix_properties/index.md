# Module `projectkoios.physkit.numerics.linear_algebra.matrix_properties`

## Current responsibility

This module owns the unit-free `matrix_is_hermitian` predicate for finite
binary64 real or complex NumPy/SciPy matrices. Dense input delegates to SciPy's
Hermiticity predicate. Sparse input compares canonical CSR structure with its
conjugate transpose without materializing a dense matrix.

The predicate requires finite nonnegative built-in-float tolerances and supports
NumPy `allclose`-style absolute and relative decisions. It owns numerical
behavior only; quantity retention and acceptance policy belong to callers.

## Local mapping

- Code: `src/python/projectkoios/physkit/numerics/linear_algebra/matrix_properties.py`
- Tests: `tests/projectkoios/physkit/numerics/linear_algebra/matrix_properties/test__matrix_is_hermitian.py`

## Navigation

- [Parent package](../index.md)
