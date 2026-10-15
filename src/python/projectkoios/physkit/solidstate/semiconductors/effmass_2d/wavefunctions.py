"""Normalized plane-wave sampling on a periodic primitive cell."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.actions import DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.model import (
    PeriodicFreeParticle2DSpectrum,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicPlaneWave2DSamplingRequest(DataObject):
    """Request normalized Bloch plane waves at fractional coordinates."""

    spectrum: PeriodicFreeParticle2DSpectrum
    fractional_coordinates: MatrixQuantity

    def __post_init__(self) -> None:
        """Check every sampling-request argument."""
        self._check_arg_spectrum()
        self._check_arg_fractional_coordinates()

    def _check_arg_spectrum(self) -> None:
        if not isinstance(self.spectrum, PeriodicFreeParticle2DSpectrum):
            raise TypeError("spectrum must be PeriodicFreeParticle2DSpectrum")

    def _check_arg_fractional_coordinates(self) -> None:
        if not isinstance(self.fractional_coordinates, MatrixQuantity):
            raise TypeError("fractional_coordinates must be MatrixQuantity")
        if not isinstance(self.fractional_coordinates.unit, Unitless):
            raise ValueError("fractional_coordinates must be dimensionless")
        if self.fractional_coordinates.magnitude.ndim != 2:
            raise ValueError("fractional_coordinates must be a matrix")
        if self.fractional_coordinates.magnitude.shape[1] != 2:
            raise ValueError("fractional_coordinates must have shape (point_count, 2)")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicPlaneWave2DSampling(ResultsObject):
    """Return physical positions and normalized mode amplitudes."""

    request: PeriodicPlaneWave2DSamplingRequest
    cartesian_positions: MatrixQuantity
    amplitudes: ComplexMatrixQuantity

    def __post_init__(self) -> None:
        """Check every sampling-result argument."""
        self._check_arg_request()
        self._check_arg_cartesian_positions()
        self._check_arg_amplitudes()

    def _check_arg_request(self) -> None:
        if not isinstance(self.request, PeriodicPlaneWave2DSamplingRequest):
            raise TypeError("request must be PeriodicPlaneWave2DSamplingRequest")

    def _check_arg_cartesian_positions(self) -> None:
        if not isinstance(self.cartesian_positions, MatrixQuantity):
            raise TypeError("cartesian_positions must be MatrixQuantity")
        expected_unit = self.request.spectrum.request.model.cell.primitive_basis.unit
        if self.cartesian_positions.unit != expected_unit:
            raise ValueError("cartesian_positions must use the model length unit")
        expected_shape = self.request.fractional_coordinates.magnitude.shape
        if self.cartesian_positions.magnitude.shape != expected_shape:
            raise ValueError("cartesian_positions must match the coordinate shape")

    def _check_arg_amplitudes(self) -> None:
        if not isinstance(self.amplitudes, ComplexMatrixQuantity):
            raise TypeError("amplitudes must be ComplexMatrixQuantity")
        point_count = self.request.fractional_coordinates.magnitude.shape[0]
        mode_count = len(self.request.spectrum.request.mode_indices)
        if self.amplitudes.magnitude.shape != (point_count, mode_count):
            raise ValueError("amplitudes shape must match points and modes")
        length_unit = (
            self.request.spectrum.request.model.cell.primitive_basis.unit.expression
        )
        expected_unit = PhysicalUnit(expression=f"1 / {length_unit}")
        if self.amplitudes.unit != expected_unit:
            raise ValueError("amplitudes must use the normalized 2D wavefunction unit")


class PeriodicPlaneWave2DSampler(
    DataObjectActionizer[
        PeriodicPlaneWave2DSamplingRequest,
        PeriodicPlaneWave2DSampling,
    ]
):
    """Sample cell-normalized Bloch plane waves in physical coordinates."""

    def action(
        self,
        *,
        request: PeriodicPlaneWave2DSamplingRequest,
    ) -> PeriodicPlaneWave2DSampling:
        """Return normalized plane-wave values at requested coordinates."""
        if not isinstance(request, PeriodicPlaneWave2DSamplingRequest):
            raise TypeError("request must be PeriodicPlaneWave2DSamplingRequest")
        spectrum = request.spectrum
        cell = spectrum.request.model.cell
        positions = (
            request.fractional_coordinates.magnitude @ cell.primitive_basis.magnitude.T
        )
        phases = positions @ spectrum.wave_vectors.magnitude.T
        normalization = 1.0 / np.sqrt(cell.area.magnitude)
        amplitudes = normalization * np.exp(1.0j * phases)
        length_unit = cell.primitive_basis.unit.expression
        return PeriodicPlaneWave2DSampling(
            request=request,
            cartesian_positions=MatrixQuantity(
                magnitude=np.asarray(positions, dtype=np.float64),
                unit=cell.primitive_basis.unit,
            ),
            amplitudes=ComplexMatrixQuantity(
                magnitude=np.asarray(amplitudes, dtype=np.complex128),
                unit=PhysicalUnit(expression=f"1 / {length_unit}"),
            ),
        )
