# `FermiDiracChemicalPotentialSolver`

## Responsibility

`FermiDiracChemicalPotentialSolver` stores deterministic bisection controls. Its
`solve(...)` façade solves an unweighted finite-level occupation target.
`solve_weighted(...)` accepts one positive unit-aware state weight per energy
and solves a target in the same state-weight unit. Both paths use the stable
`FermiDiracState` occupation evaluator; no degeneracy is implicit. The weighted
path uses particle residuals below half filling and hole residuals above half
filling to preserve precision near the ends of the inventory.

## Local mapping

- Code: `src/python/projectkoios/physkit/thermal/statmech/qm.py::FermiDiracChemicalPotentialSolver`
- Tests: `tests/projectkoios/physkit/thermal/statmech/qm/test__FermiDiracChemicalPotentialSolver__solve.py` and `test__FermiDiracChemicalPotentialSolver__solve_weighted.py`

## Navigation

- [Module](../index.md)
