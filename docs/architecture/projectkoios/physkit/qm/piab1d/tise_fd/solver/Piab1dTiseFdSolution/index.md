# `Piab1dTiseFdSolution`

## Responsibility

`Piab1dTiseFdSolution` is an immutable `ResultsObject` that retains the exact
`Piab1D`, homogeneous-Dirichlet interval, and complete unit-aware
`RealSymmetricEigendecomposition` produced by one finite-difference solve. The
decomposition owns the sparse selected-scale Hamiltonian, energy quantity, and
unitless eigenvector matrix without duplicate storage.

The compatibility accessors `hamiltonian`, `eigenpairs`, `energies`, and
`eigenvectors` return objects owned by that decomposition. Additional accessors
expose the interior-point count, spacing, and closed-form spectrum of the
represented tridiagonal matrix. `Piab1dTiseFdResults` is an exact alias for this
canonical type.

## Local mapping

- Code: `src/python/projectkoios/physkit/qm/piab1d/tise_fd/solver.py::Piab1dTiseFdSolution`
- Accessor tests: `tests/projectkoios/physkit/qm/piab1d/tise_fd/solver/test__Piab1dTiseFdSolution__energies.py`
- Spectrum tests: `tests/projectkoios/physkit/qm/piab1d/tise_fd/solver/test__Piab1dTiseFdSolution__closed_form_discrete_energies.py`
- Alias tests: `tests/projectkoios/physkit/qm/piab1d/tise_fd/solver/test__Piab1dTiseFdSolution__compatibility_alias.py`

## Navigation

- [Module](../index.md)
