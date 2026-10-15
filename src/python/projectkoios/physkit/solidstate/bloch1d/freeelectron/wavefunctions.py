"""Cell-normalized one-dimensional free-electron Bloch modes."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.actions import DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.solidstate.bloch1d.freeelectron.model import (
    FreeElectronBloch1DSpectrum,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    PhysicalUnit,
    Unitless,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class FreeElectronBloch1DWavefunctionSamplingRequest(DataObject):
    """Request cell-normalized Bloch modes at fractional coordinates."""

    spectrum: FreeElectronBloch1DSpectrum
    fractional_coordinates: VectorQuantity

    def __post_init__(self) -> None:
        """Check every sampling-request argument."""
        self._check_arg_spectrum()
        self._check_arg_fractional_coordinates()

    def _check_arg_spectrum(self) -> None:
        if not isinstance(self.spectrum, FreeElectronBloch1DSpectrum):
            raise TypeError("spectrum must be FreeElectronBloch1DSpectrum")

    def _check_arg_fractional_coordinates(self) -> None:
        if not isinstance(self.fractional_coordinates, VectorQuantity):
            raise TypeError("fractional_coordinates must be VectorQuantity")
        if not isinstance(self.fractional_coordinates.unit, Unitless):
            raise ValueError("fractional_coordinates must be dimensionless")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class FreeElectronBloch1DWavefunctionSampling(ResultsObject):
    """Return physical positions and normalized complex mode amplitudes."""

    request: FreeElectronBloch1DWavefunctionSamplingRequest
    positions: VectorQuantity
    amplitudes: ComplexMatrixQuantity

    def __post_init__(self) -> None:
        """Check every sampling-result argument."""
        self._check_arg_request()
        self._check_arg_positions()
        self._check_arg_amplitudes()

    def _check_arg_request(self) -> None:
        if not isinstance(
            self.request,
            FreeElectronBloch1DWavefunctionSamplingRequest,
        ):
            raise TypeError(
                "request must be FreeElectronBloch1DWavefunctionSamplingRequest"
            )

    def _check_arg_positions(self) -> None:
        if not isinstance(self.positions, VectorQuantity):
            raise TypeError("positions must be VectorQuantity")
        expected_unit = self.request.spectrum.request.model.cell.length.unit
        if self.positions.unit != expected_unit:
            raise ValueError("positions must use the model length unit")
        if self.positions.magnitude.shape != (
            self.request.fractional_coordinates.magnitude.size,
        ):
            raise ValueError("positions must match the requested coordinates")

    def _check_arg_amplitudes(self) -> None:
        if not isinstance(self.amplitudes, ComplexMatrixQuantity):
            raise TypeError("amplitudes must be ComplexMatrixQuantity")
        point_count = self.request.fractional_coordinates.magnitude.size
        mode_count = len(self.request.spectrum.request.mode_indices)
        if self.amplitudes.magnitude.shape != (point_count, mode_count):
            raise ValueError("amplitudes shape must match points and modes")
        length_unit = self.request.spectrum.request.model.cell.length.unit.expression
        expected_unit = PhysicalUnit(f"1 / ({length_unit} ** 0.5)")
        if self.amplitudes.unit != expected_unit:
            raise ValueError("amplitudes must use inverse-square-root length units")


class FreeElectronBloch1DWavefunctionSampler(
    DataObjectActionizer[
        FreeElectronBloch1DWavefunctionSamplingRequest,
        FreeElectronBloch1DWavefunctionSampling,
    ]
):
    """Sample exact cell-normalized free-electron Bloch modes."""

    def action(
        self,
        *,
        request: FreeElectronBloch1DWavefunctionSamplingRequest,
    ) -> FreeElectronBloch1DWavefunctionSampling:
        """Return exact plane waves at requested fractional coordinates."""
        if not isinstance(
            request,
            FreeElectronBloch1DWavefunctionSamplingRequest,
        ):
            raise TypeError(
                "request must be FreeElectronBloch1DWavefunctionSamplingRequest"
            )
        spectrum = request.spectrum
        cell = spectrum.request.model.cell
        positions = (
            request.fractional_coordinates.magnitude * cell.oriented_length.magnitude
        )
        phases = np.outer(positions, spectrum.wave_numbers.magnitude)
        normalization = 1.0 / np.sqrt(cell.length.magnitude)
        amplitudes = normalization * np.exp(1.0j * phases)
        length_unit = cell.length.unit.expression
        return FreeElectronBloch1DWavefunctionSampling(
            request=request,
            positions=VectorQuantity(
                magnitude=np.asarray(positions, dtype=np.float64),
                unit=cell.length.unit,
            ),
            amplitudes=ComplexMatrixQuantity(
                magnitude=np.asarray(amplitudes, dtype=np.complex128),
                unit=PhysicalUnit(f"1 / ({length_unit} ** 0.5)"),
            ),
        )
