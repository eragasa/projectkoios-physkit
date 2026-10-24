# Package `projectkoios.physkit.qm.piab1d.tise_fd`

## Current responsibility

This package owns the homogeneous-Dirichlet, second-order centered-difference
stationary solution for one `Piab1D` model. It retains the represented interval,
sparse Hamiltonian and complete represented real-symmetric eigensystem in
model-system units.

## Public contracts

- `Piab1dTiseFdSolver` binds one `Piab1D` at construction and accepts only the
  interior-point count when solving.
- `Piab1dTiseFdSolution` correlates that model and grid with one unit-aware
  eigensystem and exposes its Hamiltonian, energies, and eigenvectors.
- `Piab1dTiseFdResults` is an exact compatibility alias for the canonical
  solution type.

## Navigation

- [Module `solver`](solver/index.md)
- [Parent package](../index.md)
