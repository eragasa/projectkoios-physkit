# Module `projectkoios.physkit.qm.piab1d.tise_fd.solver`

## Current responsibility

This module constructs a uniform interior discretization for a solver-bound
`Piab1D`, composes its selected-scale kinetic and zero-potential operators, and
solves the complete real-symmetric finite Hamiltonian as one represented,
unit-aware eigendecomposition.

## Public classes

- [`Piab1dTiseFdSolver`](Piab1dTiseFdSolver/index.md)
- [`Piab1dTiseFdSolution`](Piab1dTiseFdSolution/index.md)

`Piab1dTiseFdResults` is an exact compatibility alias for
`Piab1dTiseFdSolution`; it does not introduce a second result type.

## Local mapping

- Code: `src/python/projectkoios/physkit/qm/piab1d/tise_fd/solver.py`
- Tests: `tests/projectkoios/physkit/qm/piab1d/tise_fd/solver/`
- Laboratory: `notebooks/qm/piab1d/piab1d.ipynb`
