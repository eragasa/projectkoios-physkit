r"""Unit-aware one-dimensional complex plane waves."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ComplexVectorQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

_METRE = PhysicalUnit(expression="meter")
_INVERSE_METRE = PhysicalUnit(expression="1 / meter")
_SECOND = PhysicalUnit(expression="second")
_INVERSE_SECOND = PhysicalUnit(expression="1 / second")


@dataclass(frozen=True, slots=True, kw_only=True)
class PlaneWave1DParameters(DataObject):
    """Represent amplitude, wave number, angular frequency, and phase offset."""

    amplitude: ScalarQuantity
    wave_number: ScalarQuantity
    angular_frequency: ScalarQuantity
    phase_offset: ScalarQuantity

    def __post_init__(self) -> None:
        """Check every one-dimensional plane-wave parameter."""
        self._check_arg_amplitude()
        self._check_arg_wave_number()
        self._check_arg_angular_frequency()
        self._check_arg_phase_offset()

    def _check_arg_amplitude(self) -> None:
        """Require a nonnegative unitless amplitude."""
        if not isinstance(self.amplitude, ScalarQuantity):
            raise TypeError("amplitude must be ScalarQuantity")
        if not isinstance(self.amplitude.unit, Unitless):
            raise ValueError("amplitude must be unitless")
        if self.amplitude.magnitude < 0.0:
            raise ValueError("amplitude must be nonnegative")

    def _check_arg_wave_number(self) -> None:
        """Require a physical inverse length."""
        if not isinstance(self.wave_number, ScalarQuantity):
            raise TypeError("wave_number must be ScalarQuantity")
        if not isinstance(self.wave_number.unit, PhysicalUnit):
            raise TypeError("wave_number must use a physical inverse-length unit")
        if not ModelSystemUnitConverter.compatible(
            self.wave_number.unit,
            _INVERSE_METRE,
        ):
            raise ValueError(
                "wave_number must use units compatible with inverse length"
            )

    def _check_arg_angular_frequency(self) -> None:
        """Require a physical inverse time."""
        if not isinstance(self.angular_frequency, ScalarQuantity):
            raise TypeError("angular_frequency must be ScalarQuantity")
        if not isinstance(self.angular_frequency.unit, PhysicalUnit):
            raise TypeError("angular_frequency must use a physical inverse-time unit")
        if not ModelSystemUnitConverter.compatible(
            self.angular_frequency.unit,
            _INVERSE_SECOND,
        ):
            raise ValueError(
                "angular_frequency must use units compatible with inverse time"
            )

    def _check_arg_phase_offset(self) -> None:
        """Require a unitless phase offset in radians."""
        if not isinstance(self.phase_offset, ScalarQuantity):
            raise TypeError("phase_offset must be ScalarQuantity")
        if not isinstance(self.phase_offset.unit, Unitless):
            raise ValueError("phase_offset must be unitless")


@dataclass(frozen=True, slots=True, kw_only=True)
class PlaneWave1DEvaluationRequest(DataObject):
    """Request one sampled one-dimensional complex plane wave."""

    model: PlaneWave1DModel
    positions: VectorQuantity
    time: ScalarQuantity

    def __post_init__(self) -> None:
        """Check every one-dimensional evaluation request argument."""
        self._check_arg_model()
        self._check_arg_positions()
        self._check_arg_time()

    def _check_arg_model(self) -> None:
        """Require a one-dimensional plane-wave model."""
        if not isinstance(self.model, PlaneWave1DModel):
            raise TypeError("model must be PlaneWave1DModel")

    def _check_arg_positions(self) -> None:
        """Require nonempty physical positions."""
        if not isinstance(self.positions, VectorQuantity):
            raise TypeError("positions must be VectorQuantity")
        if not isinstance(self.positions.unit, PhysicalUnit):
            raise TypeError("positions must use a physical length unit")
        if not ModelSystemUnitConverter.compatible(self.positions.unit, _METRE):
            raise ValueError("positions must use units compatible with length")
        if self.positions.magnitude.size == 0:
            raise ValueError("positions must be nonempty")

    def _check_arg_time(self) -> None:
        """Require one physical time."""
        if not isinstance(self.time, ScalarQuantity):
            raise TypeError("time must be ScalarQuantity")
        if not isinstance(self.time.unit, PhysicalUnit):
            raise TypeError("time must use a physical time unit")
        if not ModelSystemUnitConverter.compatible(self.time.unit, _SECOND):
            raise ValueError("time must use units compatible with time")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PlaneWave1DEvaluation(ResultsObject):
    """Retain one request and its immutable sampled complex field."""

    request: PlaneWave1DEvaluationRequest
    field_values: ComplexVectorQuantity

    def __post_init__(self) -> None:
        """Check every one-dimensional evaluation response argument."""
        self._check_arg_request()
        self._check_arg_field_values()

    def _check_arg_request(self) -> None:
        """Require the complete evaluation request."""
        if not isinstance(self.request, PlaneWave1DEvaluationRequest):
            raise TypeError("request must be PlaneWave1DEvaluationRequest")

    def _check_arg_field_values(self) -> None:
        """Require one unitless complex value per position."""
        if not isinstance(self.field_values, ComplexVectorQuantity):
            raise TypeError("field_values must be ComplexVectorQuantity")
        if not isinstance(self.field_values.unit, Unitless):
            raise ValueError("field_values must be unitless")
        if self.field_values.magnitude.shape != self.request.positions.magnitude.shape:
            raise ValueError("field_values must match positions")


@dataclass(frozen=True, slots=True)
class PlaneWave1DEvaluator:
    r"""Evaluate $A\exp[i(kx-\omega t+\phi)]$."""

    def action(self, *, request: PlaneWave1DEvaluationRequest) -> PlaneWave1DEvaluation:
        """Evaluate the phase after compatible space and time conversion."""
        if not isinstance(request, PlaneWave1DEvaluationRequest):
            raise TypeError("request must be PlaneWave1DEvaluationRequest")
        parameters = request.model.parameters
        position_unit = PhysicalUnit(
            expression=f"1 / ({parameters.wave_number.unit.expression})"
        )
        time_unit = PhysicalUnit(
            expression=f"1 / ({parameters.angular_frequency.unit.expression})"
        )
        positions = ModelSystemUnitConverter.convert_vector(
            request.positions,
            position_unit,
        )
        time = ModelSystemUnitConverter.convert_scalar(request.time, time_unit)
        phase = (
            parameters.wave_number.magnitude * positions.magnitude
            - parameters.angular_frequency.magnitude * time.magnitude
            + parameters.phase_offset.magnitude
        )
        return PlaneWave1DEvaluation(
            request=request,
            field_values=ComplexVectorQuantity(
                magnitude=np.asarray(
                    parameters.amplitude.magnitude * np.exp(1.0j * phase),
                    dtype=np.complex128,
                ),
                unit=Unitless(),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PlaneWave1DModel(
    ActionizedDataObject[PlaneWave1DEvaluationRequest, PlaneWave1DEvaluation]
):
    """Represent a unitless complex plane wave on a physical line."""

    parameters: PlaneWave1DParameters
    actionizer: PlaneWave1DEvaluator = field(
        default_factory=PlaneWave1DEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the complete parameter record."""
        self._check_arg_parameters()

    def _check_arg_parameters(self) -> None:
        """Require one-dimensional plane-wave parameters."""
        if not isinstance(self.parameters, PlaneWave1DParameters):
            raise TypeError("parameters must be PlaneWave1DParameters")

    def evaluate(
        self,
        *,
        positions: VectorQuantity,
        time: ScalarQuantity,
    ) -> PlaneWave1DEvaluation:
        """Evaluate the complex field at explicit positions and time."""
        return self._respond(
            request=PlaneWave1DEvaluationRequest(
                model=self,
                positions=positions,
                time=time,
            )
        )
