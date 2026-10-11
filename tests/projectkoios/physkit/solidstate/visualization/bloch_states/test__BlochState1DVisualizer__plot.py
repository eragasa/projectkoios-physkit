"""Tests for one-dimensional sampled Bloch-state visualization."""

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from projectkoios.physkit.solidstate.visualization.bloch_states import (
    BlochState1DPlotData,
    BlochState1DVisualizer,
)
from projectkoios.physkit.units import (
    ComplexVectorQuantity,
    PhysicalUnit,
    ScalarQuantity,
    VectorQuantity,
)


def _plane_wave_data() -> BlochState1DPlotData:
    """Return a normalized nonzero-k plane wave on one endpoint-excluded cell."""
    sample_count = 32
    length = 4.0
    positions = np.linspace(0.0, length, sample_count, endpoint=False)
    wave_number = 0.37 * np.pi / length
    amplitudes = np.exp(1j * wave_number * positions) / np.sqrt(length)
    return BlochState1DPlotData(
        positions=VectorQuantity(
            magnitude=positions,
            unit=PhysicalUnit("angstrom"),
        ),
        amplitudes=ComplexVectorQuantity(
            magnitude=amplitudes,
            unit=PhysicalUnit("1 / angstrom ** 0.5"),
        ),
        cell_length=ScalarQuantity(
            magnitude=length,
            unit=PhysicalUnit("angstrom"),
        ),
        bloch_wave_number=ScalarQuantity(
            magnitude=float(wave_number),
            unit=PhysicalUnit("1 / angstrom"),
        ),
    )


def test__plot_state__separates_amplitude_density_and_phase_units() -> None:
    """Amplitude and probability density occupy distinct labeled axes."""
    data = _plane_wave_data()
    figure, created_axes = plt.subplots(1, 3)
    axes = (created_axes[0], created_axes[1], created_axes[2])

    returned = BlochState1DVisualizer(data=data).plot_state(axes=axes)

    assert returned == axes
    assert len(axes[0].lines) == 2
    assert "1 / angstrom ** 0.5" in axes[0].get_ylabel()
    assert "** 2" in axes[1].get_ylabel()
    expected_phase = np.unwrap(np.angle(data.amplitudes.magnitude))
    np.testing.assert_allclose(axes[2].lines[0].get_ydata(), expected_phase)
    plt.close(figure)


def test__plot_cell_periodic_fourier_weights__removes_bloch_carrier() -> None:
    """A full-state plane wave becomes the zero cell-periodic Fourier mode."""
    data = _plane_wave_data()
    weights = data.cell_periodic_fourier_weights.magnitude
    indices = data.reciprocal_indices.magnitude
    figure, axis = plt.subplots()

    returned = BlochState1DVisualizer(
        data=data
    ).plot_cell_periodic_fourier_weights(ax=axis)

    assert returned is axis
    assert np.isclose(np.sum(weights), 1.0)
    assert indices[int(np.argmax(weights))] == 0.0
    assert np.max(weights) > 1.0 - 1.0e-12
    plt.close(figure)


def test__BlochState1DPlotData__rejects_a_duplicated_periodic_endpoint() -> None:
    """One period must contain unique endpoint-excluded representatives."""
    length = 4.0
    positions = np.linspace(0.0, length, 8, endpoint=True)

    with pytest.raises(ValueError, match="periodic endpoint excluded"):
        BlochState1DPlotData(
            positions=VectorQuantity(
                magnitude=positions,
                unit=PhysicalUnit("angstrom"),
            ),
            amplitudes=ComplexVectorQuantity(
                magnitude=np.ones(8, dtype=np.complex128) / np.sqrt(length),
                unit=PhysicalUnit("1 / angstrom ** 0.5"),
            ),
            cell_length=ScalarQuantity(
                magnitude=length,
                unit=PhysicalUnit("angstrom"),
            ),
            bloch_wave_number=ScalarQuantity(
                magnitude=0.0,
                unit=PhysicalUnit("1 / angstrom"),
            ),
        )


def test__BlochState1DPlotData__rejects_nonuniform_positions() -> None:
    """Discrete Fourier labels require a uniform endpoint-excluded grid."""
    with pytest.raises(ValueError, match="positions must be uniformly spaced"):
        BlochState1DPlotData(
            positions=VectorQuantity(
                magnitude=np.array([0.0, 0.2, 0.7]),
                unit=PhysicalUnit("angstrom"),
            ),
            amplitudes=ComplexVectorQuantity(
                magnitude=np.ones(3, dtype=np.complex128),
                unit=PhysicalUnit("1 / angstrom ** 0.5"),
            ),
            cell_length=ScalarQuantity(
                magnitude=0.9,
                unit=PhysicalUnit("angstrom"),
            ),
            bloch_wave_number=ScalarQuantity(
                magnitude=0.0,
                unit=PhysicalUnit("1 / angstrom"),
            ),
        )
