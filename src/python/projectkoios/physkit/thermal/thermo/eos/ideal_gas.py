r"""Unit-aware ideal-gas equation of state."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.constants import SI
from projectkoios.physkit.core.actions import ActionizedDataObject
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    VectorQuantity,
)

_KELVIN = PhysicalUnit(expression="kelvin")
_CUBIC_METRE = PhysicalUnit(expression="meter ** 3")
_MOLE = PhysicalUnit(expression="mole")
_PASCAL = PhysicalUnit(expression="pascal")


@dataclass(frozen=True, slots=True, kw_only=True)
class IdealGasPressureEvaluationRequest(DataObject):
    """Request pressures for fixed temperature over represented volumes."""

    model: IdealGasEquationOfState
    volumes: VectorQuantity
    temperature: ScalarQuantity
    pressure_unit: PhysicalUnit

    def __post_init__(self) -> None:
        """Check each ideal-gas pressure request argument."""
        self._check_arg_model()
        self._check_arg_volumes()
        self._check_arg_temperature()
        self._check_arg_pressure_unit()

    def _check_arg_model(self) -> None:
        """Require the ideal-gas equation-of-state model."""
        if not isinstance(self.model, IdealGasEquationOfState):
            raise TypeError("model must be IdealGasEquationOfState")

    def _check_arg_volumes(self) -> None:
        """Require nonempty positive physical volumes."""
        if not isinstance(self.volumes, VectorQuantity):
            raise TypeError("volumes must be VectorQuantity")
        if not isinstance(self.volumes.unit, PhysicalUnit):
            raise TypeError("volumes must use a physical volume unit")
        if not ModelSystemUnitConverter.compatible(
            self.volumes.unit,
            _CUBIC_METRE,
        ):
            raise ValueError("volumes must use units compatible with volume")
        if self.volumes.magnitude.size == 0:
            raise ValueError("volumes must be nonempty")
        if np.any(self.volumes.magnitude <= 0.0):
            raise ValueError("volumes must be positive")

    def _check_arg_temperature(self) -> None:
        """Require a physical temperature above absolute zero."""
        if not isinstance(self.temperature, ScalarQuantity):
            raise TypeError("temperature must be ScalarQuantity")
        if not isinstance(self.temperature.unit, PhysicalUnit):
            raise TypeError("temperature must use a physical temperature unit")
        if not ModelSystemUnitConverter.compatible(
            self.temperature.unit,
            _KELVIN,
        ):
            raise ValueError("temperature must use units compatible with temperature")
        if self.temperature_in_kelvin.magnitude <= 0.0:
            raise ValueError("temperature must be above absolute zero")

    def _check_arg_pressure_unit(self) -> None:
        """Require an explicit physical output pressure unit."""
        if not isinstance(self.pressure_unit, PhysicalUnit):
            raise TypeError("pressure_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.pressure_unit,
            _PASCAL,
        ):
            raise ValueError("pressure_unit must be compatible with pressure")

    @property
    def volumes_in_cubic_metres(self) -> VectorQuantity:
        """Return represented volumes in cubic metres."""
        return ModelSystemUnitConverter.convert_vector(
            quantity=self.volumes,
            target=_CUBIC_METRE,
        )

    @property
    def temperature_in_kelvin(self) -> ScalarQuantity:
        """Return absolute temperature in kelvin, including affine conversion."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.temperature,
            target=_KELVIN,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class IdealGasPressureEvaluation(ResultsObject):
    """Correlate an ideal-gas request with represented pressures."""

    request: IdealGasPressureEvaluationRequest
    pressures: VectorQuantity

    def __post_init__(self) -> None:
        """Check each ideal-gas pressure response argument."""
        self._check_arg_request()
        self._check_arg_pressures()

    def _check_arg_request(self) -> None:
        """Require the complete ideal-gas pressure request."""
        if not isinstance(self.request, IdealGasPressureEvaluationRequest):
            raise TypeError("request must be IdealGasPressureEvaluationRequest")

    def _check_arg_pressures(self) -> None:
        """Require one positive pressure per represented volume."""
        if not isinstance(self.pressures, VectorQuantity):
            raise TypeError("pressures must be VectorQuantity")
        if self.pressures.unit != self.request.pressure_unit:
            raise ValueError("pressures must use the requested pressure unit")
        if self.pressures.magnitude.shape != self.request.volumes.magnitude.shape:
            raise ValueError("pressures must match volumes")
        if np.any(self.pressures.magnitude <= 0.0):
            raise ValueError("pressures must be positive")

    @property
    def volumes(self) -> VectorQuantity:
        """Return the represented input volumes."""
        return self.request.volumes

    @property
    def temperature(self) -> ScalarQuantity:
        """Return the represented fixed temperature."""
        return self.request.temperature


@dataclass(frozen=True, slots=True)
class IdealGasPressureEvaluator:
    r"""Evaluate $P=nR_gT/V$ for one typed request."""

    def action(
        self,
        *,
        request: IdealGasPressureEvaluationRequest,
    ) -> IdealGasPressureEvaluation:
        """Return ideal-gas pressures in the requested unit."""
        if not isinstance(request, IdealGasPressureEvaluationRequest):
            raise TypeError("request must be IdealGasPressureEvaluationRequest")
        pressure_pascal = VectorQuantity(
            magnitude=(
                request.model.amount_in_moles.magnitude
                * SI.R_g
                * request.temperature_in_kelvin.magnitude
                / request.volumes_in_cubic_metres.magnitude
            ),
            unit=_PASCAL,
        )
        return IdealGasPressureEvaluation(
            request=request,
            pressures=ModelSystemUnitConverter.convert_vector(
                quantity=pressure_pascal,
                target=request.pressure_unit,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class IdealGasIsothermsEvaluationRequest(DataObject):
    """Request an ordered family of ideal-gas pressure isotherms."""

    model: IdealGasEquationOfState
    volumes: VectorQuantity
    temperatures: VectorQuantity
    pressure_unit: PhysicalUnit

    def __post_init__(self) -> None:
        """Check each ideal-gas isotherm-family request argument."""
        self._check_arg_model()
        self._check_arg_volumes()
        self._check_arg_temperatures()
        self._check_arg_pressure_unit()

    def _check_arg_model(self) -> None:
        """Require the ideal-gas equation-of-state model."""
        if not isinstance(self.model, IdealGasEquationOfState):
            raise TypeError("model must be IdealGasEquationOfState")

    def _check_arg_volumes(self) -> None:
        """Require nonempty positive physical volumes."""
        if not isinstance(self.volumes, VectorQuantity):
            raise TypeError("volumes must be VectorQuantity")
        if not isinstance(self.volumes.unit, PhysicalUnit):
            raise TypeError("volumes must use a physical volume unit")
        if not ModelSystemUnitConverter.compatible(
            self.volumes.unit,
            _CUBIC_METRE,
        ):
            raise ValueError("volumes must use units compatible with volume")
        if self.volumes.magnitude.size == 0:
            raise ValueError("volumes must be nonempty")
        if np.any(self.volumes.magnitude <= 0.0):
            raise ValueError("volumes must be positive")

    def _check_arg_temperatures(self) -> None:
        """Require nonempty physical temperatures above absolute zero."""
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
        converted = ModelSystemUnitConverter.convert_vector(
            quantity=self.temperatures,
            target=_KELVIN,
        )
        if np.any(converted.magnitude <= 0.0):
            raise ValueError("temperatures must be above absolute zero")

    def _check_arg_pressure_unit(self) -> None:
        """Require an explicit physical output pressure unit."""
        if not isinstance(self.pressure_unit, PhysicalUnit):
            raise TypeError("pressure_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.pressure_unit,
            _PASCAL,
        ):
            raise ValueError("pressure_unit must be compatible with pressure")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class IdealGasIsotherms(ResultsObject):
    """Retain an ordered iterable of ideal-gas pressure isotherms."""

    request: IdealGasIsothermsEvaluationRequest
    isotherms: tuple[IdealGasPressureEvaluation, ...]

    def __post_init__(self) -> None:
        """Check each ideal-gas isotherm-family response argument."""
        self._check_arg_request()
        self._check_arg_isotherms()

    def _check_arg_request(self) -> None:
        """Require the complete isotherm-family request."""
        if not isinstance(self.request, IdealGasIsothermsEvaluationRequest):
            raise TypeError("request must be IdealGasIsothermsEvaluationRequest")

    def _check_arg_isotherms(self) -> None:
        """Require one correlated isotherm per requested temperature."""
        if not isinstance(self.isotherms, tuple):
            raise TypeError("isotherms must be a tuple")
        if len(self.isotherms) != self.request.temperatures.magnitude.size:
            raise ValueError("isotherms must match temperatures")
        if not all(
            isinstance(isotherm, IdealGasPressureEvaluation)
            for isotherm in self.isotherms
        ):
            raise TypeError("isotherms must contain IdealGasPressureEvaluation")
        if any(
            isotherm.request.model is not self.request.model
            for isotherm in self.isotherms
        ):
            raise ValueError("isotherms must retain the requested model")
        if any(
            isotherm.volumes is not self.request.volumes for isotherm in self.isotherms
        ):
            raise ValueError("isotherms must retain the requested volumes")
        if any(
            isotherm.pressures.unit != self.request.pressure_unit
            for isotherm in self.isotherms
        ):
            raise ValueError("isotherms must use the requested pressure unit")
        represented_temperatures = np.array(
            [isotherm.temperature.magnitude for isotherm in self.isotherms],
            dtype=np.float64,
        )
        if not np.array_equal(
            represented_temperatures,
            self.request.temperatures.magnitude,
        ):
            raise ValueError("isotherms must retain requested temperature order")
        if any(
            isotherm.temperature.unit != self.request.temperatures.unit
            for isotherm in self.isotherms
        ):
            raise ValueError("isotherms must retain requested temperature units")

    def __iter__(self) -> Iterator[IdealGasPressureEvaluation]:
        """Iterate over isotherms in requested-temperature order."""
        return iter(self.isotherms)

    def __len__(self) -> int:
        """Return the number of represented temperatures."""
        return len(self.isotherms)


@dataclass(frozen=True, slots=True)
class IdealGasIsothermsEvaluator:
    """Evaluate an ordered family of ideal-gas pressure isotherms."""

    def action(
        self,
        *,
        request: IdealGasIsothermsEvaluationRequest,
    ) -> IdealGasIsotherms:
        """Return one pressure isotherm per requested temperature."""
        if not isinstance(request, IdealGasIsothermsEvaluationRequest):
            raise TypeError("request must be IdealGasIsothermsEvaluationRequest")
        return IdealGasIsotherms(
            request=request,
            isotherms=tuple(
                request.model.evaluate_pressure(
                    volumes=request.volumes,
                    temperature=ScalarQuantity(
                        magnitude=float(temperature),
                        unit=request.temperatures.unit,
                    ),
                    pressure_unit=request.pressure_unit,
                )
                for temperature in request.temperatures.magnitude
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class IdealGasEquationOfState(
    ActionizedDataObject[
        IdealGasPressureEvaluationRequest,
        IdealGasPressureEvaluation,
    ]
):
    """Represent an ideal gas with an explicit amount of substance."""

    amount: ScalarQuantity
    actionizer: IdealGasPressureEvaluator = field(
        default_factory=IdealGasPressureEvaluator,
        init=False,
        repr=False,
        compare=False,
    )
    isotherm_actionizer: IdealGasIsothermsEvaluator = field(
        default_factory=IdealGasIsothermsEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the amount-of-substance argument."""
        self._check_arg_amount()

    def _check_arg_amount(self) -> None:
        """Require a positive physical amount of substance."""
        if not isinstance(self.amount, ScalarQuantity):
            raise TypeError("amount must be ScalarQuantity")
        if not isinstance(self.amount.unit, PhysicalUnit):
            raise TypeError("amount must use a physical amount unit")
        if not ModelSystemUnitConverter.compatible(self.amount.unit, _MOLE):
            raise ValueError(
                "amount must use units compatible with amount of substance"
            )
        if self.amount_in_moles.magnitude <= 0.0:
            raise ValueError("amount must be positive")

    @property
    def amount_in_moles(self) -> ScalarQuantity:
        """Return amount of substance in moles."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.amount,
            target=_MOLE,
        )

    def evaluate_pressure(
        self,
        *,
        volumes: VectorQuantity,
        temperature: ScalarQuantity,
        pressure_unit: PhysicalUnit,
    ) -> IdealGasPressureEvaluation:
        """Evaluate one isotherm over explicit physical volumes."""
        return self._respond(
            request=IdealGasPressureEvaluationRequest(
                model=self,
                volumes=volumes,
                temperature=temperature,
                pressure_unit=pressure_unit,
            )
        )

    def evaluate_isotherms(
        self,
        *,
        volumes: VectorQuantity,
        temperatures: VectorQuantity,
        pressure_unit: PhysicalUnit,
    ) -> IdealGasIsotherms:
        """Evaluate an ordered family of fixed-temperature pressure curves."""
        return self.isotherm_actionizer.action(
            request=IdealGasIsothermsEvaluationRequest(
                model=self,
                volumes=volumes,
                temperatures=temperatures,
                pressure_unit=pressure_unit,
            )
        )
