# `MatrixVerifier`

## Responsibility

`MatrixVerifier` retains one immutable finite dense or CSR real/complex matrix
representation and evaluates generic matrix properties. Existing PhysKit
quantity inputs are retained by identity. Raw inputs must use exactly
`float64` or `complex128` storage and become canonical unitless quantities.

`is_hermitian(rtol=0.0, atol=0.0)` requires finite nonnegative built-in floats.
It compares sparse matrices with their conjugate transposes without dense
materialization.

## Local mapping

- Code: `src/python/projectkoios/physkit/linalg/matrices.py::MatrixVerifier`
- Constructor tests: `tests/projectkoios/physkit/linalg/matrices/test__MatrixVerifier__init__.py`
- Hermiticity tests: `tests/projectkoios/physkit/linalg/matrices/test__MatrixVerifier__is_hermitian.py`

## Navigation

- [Module](../index.md)
