# `BlochState1DVisualizer`

## Responsibility

`BlochState1DVisualizer` renders real and imaginary amplitudes, probability
density, unwrapped phase, and normalized cell-periodic Fourier weights. It keeps
amplitude and density on separate axes because their units differ.

The visualizer accepts or creates Matplotlib axes, returns them, and leaves
figure lifecycle and layout to the caller. It does not call `show` or introduce
an `ipywidgets` dependency.

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/visualization/bloch_states.py::BlochState1DVisualizer`
- Tests: `tests/projectkoios/physkit/solidstate/visualization/bloch_states/test__BlochState1DVisualizer__plot.py`

## Navigation

- [Module](../index.md)
