r"""Unit-aware two-dimensional complex plane waves."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
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


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PlaneWave2DParameters(DataObject):
    """Represent amplitude, wave vector, angular frequency, and phase offset."""

    amplitude: ScalarQuantity
    wave_vector: VectorQuantity
    angular_frequency: ScalarQuantity
    phase_offset: ScalarQuantity

    def __post_init__(self) -> None:
        """Check every two-dimensional plane-wave parameter."""
        self._check_arg_amplitude()
        self._check_arg_wave_vector()
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

    def _check_arg_wave_vector(self) -> None:
        """Require exactly two physical inverse-length components."""
        if not isinstance(self.wave_vector, VectorQuantity):
            raise TypeError("wave_vector must be VectorQuantity")
        if not isinstance(self.wave_vector.unit, PhysicalUnit):
            raise TypeError("wave_vector must use a physical inverse-length unit")
        if not ModelSystemUnitConverter.compatible(
            self.wave_vector.unit,
            _INVERSE_METRE,
        ):
            raise ValueError(
                "wave_vector must use units compatible with inverse length"
            )
        if self.wave_vector.magnitude.shape != (2,):
            raise ValueError("wave_vector must contain two Cartesian components")

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
class PlaneWave2DEvaluationRequest(DataObject):
    """Request one sampled two-dimensional complex plane wave."""

    model: PlaneWave2DModel
    x_coordinates: VectorQuantity
    y_coordinates: VectorQuantity
    time: ScalarQuantity

    def __post_init__(self) -> None:
        """Check every two-dimensional evaluation request argument."""
        self._check_arg_model()
        self._check_arg_x_coordinates()
        self._check_arg_y_coordinates()
        self._check_arg_time()

    def _check_arg_model(self) -> None:
        """Require a two-dimensional plane-wave model."""
        if not isinstance(self.model, PlaneWave2DModel):
            raise TypeError("model must be PlaneWave2DModel")

    def _check_arg_x_coordinates(self) -> None:
        """Require nonempty physical x coordinates."""
        self._check_coordinates(
            coordinates=self.x_coordinates,
            argument_name="x_coordinates",
        )

    def _check_arg_y_coordinates(self) -> None:
        """Require nonempty physical y coordinates."""
        self._check_coordinates(
            coordinates=self.y_coordinates,
            argument_name="y_coordinates",
        )

    @staticmethod
    def _check_coordinates(
        *,
        coordinates: VectorQuantity,
        argument_name: str,
    ) -> None:
        """Apply repeated mechanical coordinate checks."""
        if not isinstance(coordinates, VectorQuantity):
            raise TypeError(f"{argument_name} must be VectorQuantity")
        if not isinstance(coordinates.unit, PhysicalUnit):
            raise TypeError(f"{argument_name} must use a physical length unit")
        if not ModelSystemUnitConverter.compatible(coordinates.unit, _METRE):
            raise ValueError(f"{argument_name} must use units compatible with length")
        if coordinates.magnitude.size == 0:
            raise ValueError(f"{argument_name} must be nonempty")

    def _check_arg_time(self) -> None:
        """Require one physical time."""
        if not isinstance(self.time, ScalarQuantity):
            raise TypeError("time must be ScalarQuantity")
        if not isinstance(self.time.unit, PhysicalUnit):
            raise TypeError("time must use a physical time unit")
        if not ModelSystemUnitConverter.compatible(self.time.unit, _SECOND):
            raise ValueError("time must use units compatible with time")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PlaneWave2DEvaluation(ResultsObject):
    """Retain one request and its immutable sampled complex field."""

    request: PlaneWave2DEvaluationRequest
    field_values: ComplexMatrixQuantity

    def __post_init__(self) -> None:
        """Check every two-dimensional evaluation response argument."""
        self._check_arg_request()
        self._check_arg_field_values()

    def _check_arg_request(self) -> None:
        """Require the complete evaluation request."""
        if not isinstance(self.request, PlaneWave2DEvaluationRequest):
            raise TypeError("request must be PlaneWave2DEvaluationRequest")

    def _check_arg_field_values(self) -> None:
        """Require one unitless complex value per Cartesian grid point."""
        if not isinstance(self.field_values, ComplexMatrixQuantity):
            raise TypeError("field_values must be ComplexMatrixQuantity")
        if not isinstance(self.field_values.unit, Unitless):
            raise ValueError("field_values must be unitless")
        expected_shape = (
            self.request.y_coordinates.magnitude.size,
            self.request.x_coordinates.magnitude.size,
        )
        if self.field_values.magnitude.shape != expected_shape:
            raise ValueError("field_values must match the Cartesian grid")


@dataclass(frozen=True, slots=True)
class PlaneWave2DEvaluator:
    r"""Evaluate $A\exp[i(\mathbf k\cdot\mathbf r-\omega t+\phi)]$."""

    def action(self, *, request: PlaneWave2DEvaluationRequest) -> PlaneWave2DEvaluation:
        """Evaluate the phase on a Cartesian product grid."""
        if not isinstance(request, PlaneWave2DEvaluationRequest):
            raise TypeError("request must be PlaneWave2DEvaluationRequest")
        parameters = request.model.parameters
        position_unit = PhysicalUnit(
            expression=f"1 / ({parameters.wave_vector.unit.expression})"
        )
        time_unit = PhysicalUnit(
            expression=f"1 / ({parameters.angular_frequency.unit.expression})"
        )
        x_coordinates = ModelSystemUnitConverter.convert_vector(
            request.x_coordinates,
            position_unit,
        )
        y_coordinates = ModelSystemUnitConverter.convert_vector(
            request.y_coordinates,
            position_unit,
        )
        time = ModelSystemUnitConverter.convert_scalar(request.time, time_unit)
        x_grid, y_grid = np.meshgrid(
            x_coordinates.magnitude,
            y_coordinates.magnitude,
            indexing="xy",
        )
        phase = (
            parameters.wave_vector.magnitude[0] * x_grid
            + parameters.wave_vector.magnitude[1] * y_grid
            - parameters.angular_frequency.magnitude * time.magnitude
            + parameters.phase_offset.magnitude
        )
        return PlaneWave2DEvaluation(
            request=request,
            field_values=ComplexMatrixQuantity(
                magnitude=np.asarray(
                    parameters.amplitude.magnitude * np.exp(1.0j * phase),
                    dtype=np.complex128,
                ),
                unit=Unitless(),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PlaneWave2DModel(
    ActionizedDataObject[PlaneWave2DEvaluationRequest, PlaneWave2DEvaluation]
):
    """Represent a unitless complex plane wave on a physical Cartesian plane."""

    parameters: PlaneWave2DParameters
    actionizer: PlaneWave2DEvaluator = field(
        default_factory=PlaneWave2DEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the complete parameter record."""
        self._check_arg_parameters()

    def _check_arg_parameters(self) -> None:
        """Require two-dimensional plane-wave parameters."""
        if not isinstance(self.parameters, PlaneWave2DParameters):
            raise TypeError("parameters must be PlaneWave2DParameters")

    def evaluate(
        self,
        *,
        x_coordinates: VectorQuantity,
        y_coordinates: VectorQuantity,
        time: ScalarQuantity,
    ) -> PlaneWave2DEvaluation:
        """Evaluate the complex field on a Cartesian product grid."""
        return self._respond(
            request=PlaneWave2DEvaluationRequest(
                model=self,
                x_coordinates=x_coordinates,
                y_coordinates=y_coordinates,
                time=time,
            )
        )
