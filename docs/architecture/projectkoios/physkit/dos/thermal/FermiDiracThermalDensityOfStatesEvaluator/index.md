# `FermiDiracThermalDensityOfStatesEvaluator`

## Responsibility

The evaluator applies the maintained stable `FermiDiracState` occupations to a
Gaussian or histogram DOS. Its `evaluate(...)` façade constructs a checked
request and returns occupation-resolved and thermal-window quantities.

Gaussian integrals use the displayed energy grid. Histogram integrals use
midpoint occupations and explicit bin widths. Neither replaces exact
sampled-state chemical-potential inversion.

## Local mapping

- Code: `src/python/projectkoios/physkit/dos/thermal.py::FermiDiracThermalDensityOfStatesEvaluator`
- Tests: `tests/projectkoios/physkit/dos/thermal/test__FermiDiracThermalDensityOfStatesEvaluator__evaluate.py`

## Navigation

- [Module](../index.md)
