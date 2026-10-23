# `Piab1dTiseFdSolver`

## Responsibility

`Piab1dTiseFdSolver(model=...)` binds one nominal `Piab1D` physical model.
`solve(interior_points=...)` then selects only the spatial discretization and
returns a `Piab1dTiseFdSolution`. The point count must be a positive built-in
`int`; booleans and integer scalar types owned by third-party libraries are not
accepted.

The construction uses a uniform grid with `interior_points + 2` total points,
homogeneous Dirichlet boundaries, the centered three-point Laplacian, zero
interior potential, and a complete unit-aware real-symmetric eigendecomposition
composed over the unit-free backend. Binding the model does not change those
represented numerical or physical choices.

## Local mapping

- Code: `src/python/projectkoios/physkit/qm/piab1d/tise_fd/solver.py::Piab1dTiseFdSolver`
- Constructor tests: `tests/projectkoios/physkit/qm/piab1d/tise_fd/solver/test__Piab1dTiseFdSolver__init__.py`
- Solve tests: `tests/projectkoios/physkit/qm/piab1d/tise_fd/solver/test__Piab1dTiseFdSolver__solve.py`

## Navigation

- [Module](../index.md)
