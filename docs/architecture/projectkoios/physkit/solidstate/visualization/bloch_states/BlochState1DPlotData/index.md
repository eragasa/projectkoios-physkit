# `BlochState1DPlotData`

## Responsibility

`BlochState1DPlotData` correlates a positive cell length, a uniform
endpoint-excluded coordinate grid spanning that period, complex
one-dimensional Bloch-state amplitudes, and an explicit Bloch wave number. It
derives probability density, the cell-periodic factor, centered
integer reciprocal indices, and normalized Fourier weights.

Physical data require length coordinates, inverse-square-root-length
amplitudes, and reciprocal-length wave number. Fully nondimensional data use
`Unitless` throughout. The record validates but never renormalizes amplitudes.

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/visualization/bloch_states.py::BlochState1DPlotData`
- Tests: `tests/projectkoios/physkit/solidstate/visualization/bloch_states/test__BlochState1DVisualizer__plot.py`

## Navigation

- [Module](../index.md)
