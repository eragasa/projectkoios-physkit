# Module `projectkoios.physkit.thermal.statmech.qm`

## Current responsibility

This module owns reusable finite-level canonical and Fermi–Dirac occupation
analysis, explicitly weighted chemical-potential inversion, and
occupation-weighted sampled-eigenfunction density evaluation. Energy,
thermal-energy, chemical-potential, and state-weight units are explicit.
Degeneracy is never introduced implicitly.

## Supported façades and responses

- [`FermiDiracState`](FermiDiracState/index.md)
- [`FermiDiracOccupationEvaluation`](FermiDiracOccupationEvaluation/index.md)
- [`FermiDiracChemicalPotentialSolver`](FermiDiracChemicalPotentialSolver/index.md)
- [`FermiDiracChemicalPotentialResult`](FermiDiracChemicalPotentialResult/index.md)
- [`WeightedFermiDiracChemicalPotentialResult`](WeightedFermiDiracChemicalPotentialResult/index.md)
- [`CanonicalLevelProbabilityModel`](CanonicalLevelProbabilityModel/index.md)
- [`CanonicalLevelProbabilityEvaluation`](CanonicalLevelProbabilityEvaluation/index.md)
- [`OccupiedEigenfunctionDensityEvaluator`](OccupiedEigenfunctionDensityEvaluator/index.md)
- [`OccupiedEigenfunctionDensityEvaluation`](OccupiedEigenfunctionDensityEvaluation/index.md)

## Typed operation requests

- [`FermiDiracOccupationEvaluationRequest`](FermiDiracOccupationEvaluationRequest/index.md)
- [`FermiDiracChemicalPotentialSolveRequest`](FermiDiracChemicalPotentialSolveRequest/index.md)
- [`WeightedFermiDiracChemicalPotentialSolveRequest`](WeightedFermiDiracChemicalPotentialSolveRequest/index.md)
- [`CanonicalLevelProbabilityEvaluationRequest`](CanonicalLevelProbabilityEvaluationRequest/index.md)
- [`OccupiedEigenfunctionDensityEvaluationRequest`](OccupiedEigenfunctionDensityEvaluationRequest/index.md)

Internal ActionObjects implement the typed `action(request=...)` boundary.
Ordinary callers use façade methods such as `evaluate_occupations(...)`,
`solve(...)`, and `evaluate(...)`.

## Local mapping

- Code: `src/python/projectkoios/physkit/thermal/statmech/qm.py`
- Tests: `tests/projectkoios/physkit/thermal/statmech/qm/`
- Laboratory: `notebooks/qm/piab1d/piab1d__fermions.ipynb`
