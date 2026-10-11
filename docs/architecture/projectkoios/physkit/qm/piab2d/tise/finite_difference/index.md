# Package `projectkoios.physkit.qm.piab2d.tise.finite_difference`

## Current responsibility

This package owns centered finite-difference stationary states for
`ObliquePiab2D`, a homogeneous-Dirichlet free particle in a direct-basis
parallelogram. It owns the full-grid request, interior unknown ordering, sparse
Hamiltonian, selected low eigensystem, constant physical area weight, and mapped
interior coordinates.

It does not own periodic seams, Bloch twists, reciprocal plane waves,
rectangular analytical solutions, potentials, or anisotropic effective mass.

## Public contracts

- `ObliquePiab2DTiseFdSolveRequest` identifies the model, two full-grid point
  counts including boundaries, and selected state count.
- `ObliquePiab2DTiseFdSolver` constructs the inverse-metric pure and mixed
  stencil, applies kinetic scaling, and solves the low real-symmetric spectrum.
- `ObliquePiab2DTiseFdSolution` retains fractional and Cartesian coordinates,
  quadrature weight, energies, normalized eigenfunction columns, and canonical
  CSR Hamiltonian.

## Navigation

- [Schematic](schematics.md)
- [Numerical contract](numeric.md)
- [Scientific scope](scientific.md)
- [Parent TISE package](../index.md)
- [PIAB2D package](../../index.md)

## Local mapping

- Source:
  `src/python/projectkoios/physkit/qm/piab2d/tise/finite_difference/`
- Model geometry: `src/python/projectkoios/physkit/qm/piab2d/oblique.py`
- Tests:
  `tests/projectkoios/physkit/qm/piab2d/tise/finite_difference/`
- Laboratory:
  `notebooks/qm/piab2d/tise/numerical/oblique-dirichlet-free-particle.ipynb`
