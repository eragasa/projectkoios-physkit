r"""Unit-aware one-dimensional sampled Bloch-state visualization."""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.units import (
    ComplexVectorQuantity,
    ModelSystemUnit,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class BlochState1DPlotData(DataObject):
    r"""Represent one sampled one-dimensional Bloch state.

    Parameters
    ----------
    positions
        Strictly increasing, uniformly spaced endpoint-excluded coordinates.
    amplitudes
        Complex Bloch-state values :math:`\psi_k(x)` at each coordinate.
    cell_length
        Positive period represented by the endpoint-excluded grid.
    bloch_wave_number
        Scalar :math:`k` in the reciprocal coordinate unit. For a fully
        nondimensional representation, all three inputs must be unitless.

    Notes
    -----
    Physical amplitudes must have inverse-square-root-length dimensions. The
    data contract checks units and grid geometry but does not silently
    renormalize the state.
    """

    positions: VectorQuantity
    amplitudes: ComplexVectorQuantity
    cell_length: ScalarQuantity
    bloch_wave_number: ScalarQuantity

    def __post_init__(self) -> None:
        """Check correlated samples, uniform spacing, and dimensional units."""
        if not isinstance(self.positions, VectorQuantity):
            raise TypeError("positions must be VectorQuantity")
        if self.positions.magnitude.size < 2:
            raise ValueError("positions must contain at least two samples")
        spacings = np.diff(self.positions.magnitude)
        if np.any(spacings <= 0.0):
            raise ValueError("positions must be strictly increasing")
        spacing_scale = max(
            float(np.max(np.abs(spacings))),
            np.finfo(np.float64).tiny,
        )
        spacing_tolerance = 64.0 * np.finfo(np.float64).eps * spacing_scale
        if not np.allclose(
            spacings,
            spacings[0],
            rtol=64.0 * np.finfo(np.float64).eps,
            atol=spacing_tolerance,
        ):
            raise ValueError("positions must be uniformly spaced")
        if not isinstance(self.amplitudes, ComplexVectorQuantity):
            raise TypeError("amplitudes must be ComplexVectorQuantity")
        if self.amplitudes.magnitude.shape != self.positions.magnitude.shape:
            raise ValueError("amplitudes must match positions")
        if np.all(self.amplitudes.magnitude == 0.0):
            raise ValueError("amplitudes must contain a nonzero state")
        if not isinstance(self.cell_length, ScalarQuantity):
            raise TypeError("cell_length must be ScalarQuantity")
        if self.cell_length.magnitude <= 0.0:
            raise ValueError("cell_length must be positive")
        if not isinstance(self.bloch_wave_number, ScalarQuantity):
            raise TypeError("bloch_wave_number must be ScalarQuantity")
        self._check_units()
        self._check_endpoint_excluded_period(spacings=spacings)

    def _check_units(self) -> None:
        """Require a coherent physical or fully nondimensional convention."""
        position_unit = self.positions.unit
        amplitude_unit = self.amplitudes.unit
        wave_number_unit = self.bloch_wave_number.unit
        if isinstance(position_unit, Unitless):
            if (
                not isinstance(amplitude_unit, Unitless)
                or not isinstance(self.cell_length.unit, Unitless)
                or not isinstance(wave_number_unit, Unitless)
            ):
                raise ValueError(
                    "positions, amplitudes, cell length, and wave number must all "
                    "be unitless"
                )
            return
        if not isinstance(position_unit, PhysicalUnit):
            raise TypeError("positions must use PhysicalUnit or Unitless")
        expected_amplitude_unit = PhysicalUnit(
            f"1 / ({position_unit.expression}) ** 0.5"
        )
        if not ModelSystemUnitConverter.compatible(
            amplitude_unit,
            expected_amplitude_unit,
        ):
            raise ValueError("amplitudes must use inverse-square-root-length units")
        if not ModelSystemUnitConverter.compatible(
            self.cell_length.unit,
            position_unit,
        ):
            raise ValueError("cell_length must use length units")
        expected_wave_number_unit = PhysicalUnit(f"1 / ({position_unit.expression})")
        if not ModelSystemUnitConverter.compatible(
            wave_number_unit,
            expected_wave_number_unit,
        ):
            raise ValueError("bloch_wave_number must use reciprocal-length units")

    def _check_endpoint_excluded_period(self, *, spacings: np.ndarray) -> None:
        """Require grid spacing equal to one period divided by sample count."""
        if isinstance(self.positions.unit, Unitless):
            length = self.cell_length.magnitude
        else:
            length = ModelSystemUnitConverter.convert_scalar(
                self.cell_length,
                self.positions.unit,
            ).magnitude
        expected_spacing = length / self.positions.magnitude.size
        tolerance = (
            128.0
            * np.finfo(np.float64).eps
            * max(abs(expected_spacing), np.finfo(np.float64).tiny)
        )
        if not np.isclose(
            spacings[0],
            expected_spacing,
            rtol=128.0 * np.finfo(np.float64).eps,
            atol=tolerance,
        ):
            raise ValueError(
                "positions must sample one cell with the periodic endpoint excluded"
            )

    @property
    def probability_density(self) -> VectorQuantity:
        r"""Return :math:`|\psi_k(x)|^2` in the squared-amplitude unit."""
        if isinstance(self.amplitudes.unit, Unitless):
            density_unit: ModelSystemUnit = Unitless()
        else:
            density_unit = PhysicalUnit(f"({self.amplitudes.unit.expression}) ** 2")
        return VectorQuantity(
            magnitude=np.abs(self.amplitudes.magnitude) ** 2,
            unit=density_unit,
        )

    @property
    def cell_periodic_amplitudes(self) -> ComplexVectorQuantity:
        r"""Return :math:`u_k(x)=e^{-ikx}\psi_k(x)` on the sampled cell."""
        if isinstance(self.positions.unit, Unitless):
            wave_number = self.bloch_wave_number.magnitude
        else:
            if not isinstance(self.positions.unit, PhysicalUnit):
                raise TypeError("physical conversion requires PhysicalUnit positions")
            target_unit = PhysicalUnit(f"1 / ({self.positions.unit.expression})")
            wave_number = ModelSystemUnitConverter.convert_scalar(
                self.bloch_wave_number,
                target_unit,
            ).magnitude
        phase = np.exp(-1j * wave_number * self.positions.magnitude)
        return ComplexVectorQuantity(
            magnitude=phase * self.amplitudes.magnitude,
            unit=self.amplitudes.unit,
        )

    @property
    def reciprocal_indices(self) -> VectorQuantity:
        """Return centered integer reciprocal-mode labels for the sampled grid."""
        sample_count = self.positions.magnitude.size
        indices = np.fft.fftshift(np.fft.fftfreq(sample_count, d=1.0 / sample_count))
        return VectorQuantity(
            magnitude=np.asarray(indices, dtype=np.float64),
            unit=Unitless(),
        )

    @property
    def cell_periodic_fourier_weights(self) -> VectorQuantity:
        """Return normalized discrete Fourier weights of :math:`u_k`."""
        # Removing the Bloch carrier before the FFT makes the returned integer
        # labels reciprocal-cell modes rather than a leakage-prone transform of
        # the twisted full state.
        coefficients = np.fft.fftshift(
            np.fft.fft(self.cell_periodic_amplitudes.magnitude)
        )
        weights = np.abs(coefficients) ** 2
        normalization = float(np.sum(weights))
        return VectorQuantity(
            magnitude=np.asarray(weights / normalization, dtype=np.float64),
            unit=Unitless(),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class BlochState1DVisualizer:
    """Render amplitude, probability, phase, and reciprocal weights.

    Parameters
    ----------
    data
        Checked sampled Bloch state and wave number.
    """

    data: BlochState1DPlotData

    def __post_init__(self) -> None:
        """Require complete sampled-state plot data."""
        if not isinstance(self.data, BlochState1DPlotData):
            raise TypeError("data must be BlochState1DPlotData")

    def plot_state(
        self,
        *,
        axes: tuple[Axes, Axes, Axes] | None = None,
        title: str = "Sampled one-dimensional Bloch state",
    ) -> tuple[Axes, Axes, Axes]:
        """Plot complex amplitude, probability density, and unwrapped phase.

        Parameters
        ----------
        axes
            Optional amplitude, probability-density, and phase axes.
        title
            Figure-level scientific description applied to the amplitude axis.

        Returns
        -------
        tuple[Axes, Axes, Axes]
            Axes containing amplitude, density, and phase respectively.
        """
        if axes is None:
            _, created_axes = plt.subplots(1, 3, figsize=(13.0, 3.8))
            axes = (created_axes[0], created_axes[1], created_axes[2])
        if (
            not isinstance(axes, tuple)
            or len(axes) != 3
            or not all(isinstance(axis, Axes) for axis in axes)
        ):
            raise TypeError("axes must be a three-Axes tuple or None")
        if not isinstance(title, str):
            raise TypeError("title must be a string")

        amplitude_axis, density_axis, phase_axis = axes
        positions = self.data.positions.magnitude
        amplitudes = self.data.amplitudes.magnitude
        density = self.data.probability_density

        # Keep amplitude and probability density on separate axes because their
        # physical units differ even for a normalized state.
        amplitude_axis.plot(positions, amplitudes.real, label=r"$\mathrm{Re}\,\psi$")
        amplitude_axis.plot(positions, amplitudes.imag, label=r"$\mathrm{Im}\,\psi$")
        amplitude_axis.set(
            xlabel=f"Position ({self.data.positions.unit.expression})",
            ylabel=f"Amplitude ({self.data.amplitudes.unit.expression})",
            title=title,
        )
        amplitude_axis.legend()

        density_axis.plot(positions, density.magnitude)
        density_axis.set(
            xlabel=f"Position ({self.data.positions.unit.expression})",
            ylabel=f"Probability density ({density.unit.expression})",
            title="Probability density",
        )

        phase_axis.plot(positions, np.unwrap(np.angle(amplitudes)))
        phase_axis.set(
            xlabel=f"Position ({self.data.positions.unit.expression})",
            ylabel="Unwrapped phase (radian)",
            title="Bloch-state phase",
        )
        for axis in axes:
            axis.grid(alpha=0.25)
        return axes

    def plot_cell_periodic_fourier_weights(
        self,
        *,
        ax: Axes | None = None,
        title: str = "Cell-periodic Fourier weights",
    ) -> Axes:
        """Plot normalized reciprocal-index weights after removing Bloch phase.

        Parameters
        ----------
        ax
            Optional Matplotlib axes.
        title
            Plot title.

        Returns
        -------
        Axes
            Axes containing normalized integer reciprocal-mode weights.
        """
        if ax is None:
            _, ax = plt.subplots(figsize=(7.2, 4.2))
        if not isinstance(ax, Axes):
            raise TypeError("ax must be matplotlib.axes.Axes or None")
        if not isinstance(title, str):
            raise TypeError("title must be a string")
        ax.stem(
            self.data.reciprocal_indices.magnitude,
            self.data.cell_periodic_fourier_weights.magnitude,
        )
        ax.set(
            xlabel="Reciprocal-cell index n",
            ylabel="Normalized Fourier weight",
            title=title,
        )
        ax.grid(alpha=0.25)
        return ax
