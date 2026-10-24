# Module `projectkoios.physkit.numerics.eigenproblems.residuals`

## Current responsibility

This module owns the unit-free `eigenpair_residual_norms` kernel. Given a finite
operator, selected eigenvalues, and column eigenvectors, it returns one binary64
Euclidean norm for each column of
`operator @ eigenvectors - eigenvectors * eigenvalues`.

Dense and SciPy sparse operators may use exactly `float64` or `complex128`
storage. Accepted sparse formats are canonicalized to CSR without dense
materialization before validation and multiplication. Eigenvalue and
eigenvector arrays require the same exact dtype boundary; no
boolean, string, object, lower-precision, or integer coercion occurs.

The kernel validates finite inputs and finite arithmetic results, requires one
to the operator dimension selected columns, and correlates all dimensions. It
defines no tolerance, pass/fail decision, represented units, or scientific
acceptance policy.

## Local mapping

- Code: `src/python/projectkoios/physkit/numerics/eigenproblems/residuals.py`
- Tests: `tests/projectkoios/physkit/numerics/eigenproblems/residuals/test__eigenpair_residual_norms.py`

## Navigation

- [Parent package](../index.md)
