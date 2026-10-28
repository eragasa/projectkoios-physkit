# Module `projectkoios.physkit.numerics.eigenproblems.complex_hermitian`

## Current responsibility

This module owns complete binary64 dense complex-Hermitian eigensolver behavior
over unit-free NumPy arrays. It accepts exactly finite `complex128` matrices,
requires a nonempty square representation and exact Hermiticity, and delegates
the complete solve to `numpy.linalg.eigh`.

`ComplexHermitianEigenpairResult` retains immutable `float64` eigenvalues and
`complex128` column eigenvectors. `ComplexHermitianEigenpairSolver` returns
eigenvalues in nondecreasing order and canonicalizes each eigenvector's global
phase using its largest-magnitude component. The convention removes arbitrary
per-column phase for an ordinary eigensystem but does not define a unique basis
inside a degenerate eigenspace.

The solver does not accept sparse matrices, quantities, tolerance policy, or
scientific metadata. The represented owner
`projectkoios.physkit.linalg.eigen` supplies an unchanged matrix magnitude and
reattaches the exact matrix unit to eigenvalues and residual evidence. Sparse
complex-Hermitian solver policy remains outside this bounded dense capability.

## Local mapping

- Code: `src/python/projectkoios/physkit/numerics/eigenproblems/complex_hermitian.py`
- Tests: `tests/projectkoios/physkit/numerics/eigenproblems/complex_hermitian/`

## Navigation

- [Parent package](../index.md)
