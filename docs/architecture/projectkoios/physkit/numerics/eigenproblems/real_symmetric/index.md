# Module `projectkoios.physkit.numerics.eigenproblems.real_symmetric`

## Current responsibility

This module owns complete and lowest-state binary64 real-symmetric eigensolver
algorithms over unit-free NumPy dense and SciPy sparse matrices.

`RealSymmetricEigenpairResult` contains immutable numerical eigenvalue and
eigenvector arrays only. `RealSymmetricEigenpairSolver` validates exact symmetry
and delegates to dense or sparse SciPy algorithms while preserving canonical
ordered, sign-normalized eigenvectors. Its raw `execute_lowest` capability
remains available to numerical users; represented `linalg.eigen` construction
uses the complete solver once and performs later rank selection from stored
evidence.

The result is a numerical backend output, not the represented scientific
operator decomposition. The latter is owned by
`projectkoios.physkit.linalg.eigen`.

## Local mapping

- Code: `src/python/projectkoios/physkit/numerics/eigenproblems/real_symmetric.py`
- Tests: `tests/projectkoios/physkit/numerics/eigenproblems/real_symmetric/`

## Navigation

- [Parent package](../index.md)
