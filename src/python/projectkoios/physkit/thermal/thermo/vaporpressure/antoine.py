"""Unit-aware Antoine vapor-pressure models and evaluation."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

_KELVIN = PhysicalUnit(expression="kelvin")
_PASCAL = PhysicalUnit(expression="pascal")


@dataclass(frozen=True, slots=True, kw_only=True)
class AntoineTemperatureRange(DataObject):
    """Represent an inclusive physical validity interval for temperature."""

    lower: ScalarQuantity
    upper: ScalarQuantity

    def __post_init__(self) -> None:
        """Check each validity-range argument."""
        self._check_arg_lower()
        self._check_arg_upper()

    def _check_arg_lower(self) -> None:
        """Require a physical temperature above absolute zero."""
        if not isinstance(self.lower, ScalarQuantity):
            raise TypeError("lower must be ScalarQuantity")
        if not isinstance(self.lower.unit, PhysicalUnit):
            raise TypeError("lower must use a physical temperature unit")
        if not ModelSystemUnitConverter.compatible(self.lower.unit, _KELVIN):
            raise ValueError("lower must use units compatible with temperature")
        if self.lower_in_kelvin.magnitude <= 0.0:
            raise ValueError("lower must be above absolute zero")

    def _check_arg_upper(self) -> None:
        """Require a compatible temperature strictly above the lower bound."""
        if not isinstance(self.upper, ScalarQuantity):
            raise TypeError("upper must be ScalarQuantity")
        if not isinstance(self.upper.unit, PhysicalUnit):
            raise TypeError("upper must use a physical temperature unit")
        if not ModelSystemUnitConverter.compatible(self.upper.unit, _KELVIN):
            raise ValueError("upper must use units compatible with temperature")
        lower_kelvin = self.lower_in_kelvin.magnitude
        upper_kelvin = self.upper_in_kelvin.magnitude
        if upper_kelvin <= lower_kelvin:
            raise ValueError("upper must be greater than lower")

    @property
    def lower_in_kelvin(self) -> ScalarQuantity:
        """Return the inclusive lower bound in kelvin."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.lower,
            target=_KELVIN,
        )

    @property
    def upper_in_kelvin(self) -> ScalarQuantity:
        """Return the inclusive upper bound in kelvin."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.upper,
            target=_KELVIN,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class AntoineVaporPressureEvaluationRequest(DataObject):
    """Request vapor pressures from one model and temperature vector."""

    model: AntoineVaporPressureModel
    temperatures: VectorQuantity
    pressure_unit: PhysicalUnit
    enforce_valid_range: bool

    def __post_init__(self) -> None:
        """Check each Antoine-evaluation request argument."""
        self._check_arg_model()
        self._check_arg_temperatures()
        self._check_arg_pressure_unit()
        self._check_arg_enforce_valid_range()

    def _check_arg_model(self) -> None:
        """Require the Antoine coefficient model to evaluate."""
        if not isinstance(self.model, AntoineVaporPressureModel):
            raise TypeError("model must be AntoineVaporPressureModel")

    def _check_arg_temperatures(self) -> None:
        """Require nonempty positive physical temperatures."""
        if not isinstance(self.temperatures, VectorQuantity):
            raise TypeError("temperatures must be VectorQuantity")
        if not isinstance(self.temperatures.unit, PhysicalUnit):
            raise TypeError("temperatures must use a physical temperature unit")
        if not ModelSystemUnitConverter.compatible(
            self.temperatures.unit,
            _KELVIN,
        ):
            raise ValueError("temperatures must use units compatible with temperature")
        if self.temperatures.magnitude.size == 0:
            raise ValueError("temperatures must be nonempty")
        if np.any(self.temperatures_in(unit=_KELVIN).magnitude <= 0.0):
            raise ValueError("temperatures must be above absolute zero")

    def temperatures_in(self, *, unit: PhysicalUnit) -> VectorQuantity:
        """Return requested temperatures converted to a compatible unit."""
        if not isinstance(unit, PhysicalUnit):
            raise TypeError("unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(unit, _KELVIN):
            raise ValueError("unit must be compatible with temperature")
        return ModelSystemUnitConverter.convert_vector(
            quantity=self.temperatures,
            target=unit,
        )

    def _check_arg_pressure_unit(self) -> None:
        """Require the requested physical pressure unit."""
        if not isinstance(self.pressure_unit, PhysicalUnit):
            raise TypeError("pressure_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.pressure_unit,
            _PASCAL,
        ):
            raise ValueError("pressure_unit must be compatible with pressure")

    def _check_arg_enforce_valid_range(self) -> None:
        """Require a built-in boolean validity-range policy."""
        if type(self.enforce_valid_range) is not bool:
            raise TypeError("enforce_valid_range must be a built-in boolean")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class AntoineVaporPressureEvaluation(ResultsObject):
    """Correlate an Antoine request with evaluated vapor pressures."""

    request: AntoineVaporPressureEvaluationRequest
    pressures: VectorQuantity

    def __post_init__(self) -> None:
        """Check each Antoine-evaluation response argument."""
        self._check_arg_request()
        self._check_arg_pressures()

    def _check_arg_request(self) -> None:
        """Require the complete Antoine-evaluation request."""
        if not isinstance(self.request, AntoineVaporPressureEvaluationRequest):
            raise TypeError("request must be AntoineVaporPressureEvaluationRequest")

    def _check_arg_pressures(self) -> None:
        """Require finite nonnegative pressures matching requested temperatures."""
        if not isinstance(self.pressures, VectorQuantity):
            raise TypeError("pressures must be VectorQuantity")
        if self.pressures.unit != self.request.pressure_unit:
            raise ValueError("pressures must use the requested pressure unit")
        if self.pressures.magnitude.shape != self.request.temperatures.magnitude.shape:
            raise ValueError("pressures must match temperatures")
        if np.any(self.pressures.magnitude < 0.0):
            raise ValueError("pressures must be nonnegative")

    @property
    def temperatures(self) -> VectorQuantity:
        """Return the represented input temperatures."""
        return self.request.temperatures


@dataclass(frozen=True, slots=True)
class AntoineVaporPressureEvaluator:
    """Evaluate a typed Antoine vapor-pressure request."""

    def action(
        self,
        *,
        request: AntoineVaporPressureEvaluationRequest,
    ) -> AntoineVaporPressureEvaluation:
        """Return vapor pressures in the explicitly requested pressure unit."""
        if not isinstance(request, AntoineVaporPressureEvaluationRequest):
            raise TypeError("request must be AntoineVaporPressureEvaluationRequest")
        temperatures_kelvin = request.temperatures_in(unit=_KELVIN).magnitude
        validity_range = request.model.valid_temperature_range
        if (
            request.enforce_valid_range
            and validity_range is not None
            and (
                np.any(temperatures_kelvin < validity_range.lower_in_kelvin.magnitude)
                or np.any(
                    temperatures_kelvin > validity_range.upper_in_kelvin.magnitude
                )
            )
        ):
            raise ValueError("temperatures lie outside the model validity range")

        coefficient_temperatures = request.temperatures_in(
            unit=request.model.coefficient_temperature_unit
        ).magnitude
        denominator = coefficient_temperatures + request.model.C.magnitude
        if np.any(denominator == 0.0):
            raise ValueError("temperatures reach the Antoine denominator singularity")
        logarithmic_pressure = (
            request.model.A.magnitude - request.model.B.magnitude / denominator
        )
        with np.errstate(over="ignore", invalid="ignore"):
            native_pressure_magnitudes = np.power(10.0, logarithmic_pressure)
        if not np.all(np.isfinite(native_pressure_magnitudes)):
            raise ValueError("evaluated pressures must be finite")
        native_pressures = VectorQuantity(
            magnitude=native_pressure_magnitudes,
            unit=request.model.coefficient_pressure_unit,
        )
        pressures = ModelSystemUnitConverter.convert_vector(
            quantity=native_pressures,
            target=request.pressure_unit,
        )
        return AntoineVaporPressureEvaluation(
            request=request,
            pressures=pressures,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class AntoineVaporPressureModel(
    ActionizedDataObject[
        AntoineVaporPressureEvaluationRequest,
        AntoineVaporPressureEvaluation,
    ]
):
    r"""Represent one unit-aware Antoine vapor-pressure coefficient set.

    The represented correlation is

    .. math::

       \log_{10} P_{\mathrm{native}}
       = A - \frac{B}{T_{\mathrm{native}} + C}.
    """

    A: ScalarQuantity
    B: ScalarQuantity
    C: ScalarQuantity
    coefficient_temperature_unit: PhysicalUnit
    coefficient_pressure_unit: PhysicalUnit
    valid_temperature_range: AntoineTemperatureRange | None
    source: str
    actionizer: AntoineVaporPressureEvaluator = field(
        default_factory=AntoineVaporPressureEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check each Antoine-model argument."""
        self._check_arg_A()
        self._check_arg_B()
        self._check_arg_C()
        self._check_arg_coefficient_temperature_unit()
        self._check_arg_coefficient_pressure_unit()
        self._check_arg_valid_temperature_range()
        self._check_arg_source()

    def _check_arg_A(self) -> None:
        """Require a unitless scalar logarithmic intercept."""
        if not isinstance(self.A, ScalarQuantity):
            raise TypeError("A must be ScalarQuantity")
        if not isinstance(self.A.unit, Unitless):
            raise ValueError("A must be unitless")

    def _check_arg_B(self) -> None:
        """Require a scalar numerator in the coefficient temperature unit."""
        if not isinstance(self.B, ScalarQuantity):
            raise TypeError("B must be ScalarQuantity")
        if self.B.unit != self.coefficient_temperature_unit:
            raise ValueError("B must use the coefficient temperature unit")

    def _check_arg_C(self) -> None:
        """Require a scalar offset in the coefficient temperature unit."""
        if not isinstance(self.C, ScalarQuantity):
            raise TypeError("C must be ScalarQuantity")
        if self.C.unit != self.coefficient_temperature_unit:
            raise ValueError("C must use the coefficient temperature unit")

    def _check_arg_coefficient_temperature_unit(self) -> None:
        """Require the physical temperature unit used by the coefficient table."""
        if not isinstance(self.coefficient_temperature_unit, PhysicalUnit):
            raise TypeError("coefficient_temperature_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.coefficient_temperature_unit,
            _KELVIN,
        ):
            raise ValueError(
                "coefficient_temperature_unit must be compatible with temperature"
            )

    def _check_arg_coefficient_pressure_unit(self) -> None:
        """Require the physical pressure unit produced by the coefficient table."""
        if not isinstance(self.coefficient_pressure_unit, PhysicalUnit):
            raise TypeError("coefficient_pressure_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.coefficient_pressure_unit,
            _PASCAL,
        ):
            raise ValueError(
                "coefficient_pressure_unit must be compatible with pressure"
            )

    def _check_arg_valid_temperature_range(self) -> None:
        """Require an optional unit-aware inclusive validity interval."""
        if self.valid_temperature_range is not None and not isinstance(
            self.valid_temperature_range,
            AntoineTemperatureRange,
        ):
            raise TypeError(
                "valid_temperature_range must be AntoineTemperatureRange or None"
            )

    def _check_arg_source(self) -> None:
        """Require a nonempty coefficient source or provenance description."""
        if type(self.source) is not str:
            raise TypeError("source must be a string")
        if not self.source:
            raise ValueError("source must be nonempty")

    def evaluate(
        self,
        *,
        temperatures: VectorQuantity,
        pressure_unit: PhysicalUnit,
        enforce_valid_range: bool = True,
    ) -> AntoineVaporPressureEvaluation:
        """Evaluate vapor pressure for explicit temperatures and output units."""
        return self._respond(
            request=AntoineVaporPressureEvaluationRequest(
                model=self,
                temperatures=temperatures,
                pressure_unit=pressure_unit,
                enforce_valid_range=enforce_valid_range,
            )
        )
