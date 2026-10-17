r"""Unit-aware thermal de Broglie wavelength and degeneracy parameter."""

from __future__ import annotations

import math
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
    Unitless,
    VectorQuantity,
)

_KELVIN = PhysicalUnit(expression="kelvin")
_KILOGRAM = PhysicalUnit(expression="kilogram")
_METRE = PhysicalUnit(expression="meter")
_NUMBER_DENSITY_SI = PhysicalUnit(expression="meter ** -3")


@dataclass(frozen=True, slots=True, kw_only=True)
class ThermalDeBroglieWavelengthEvaluationRequest(DataObject):
    """Request thermal wavelengths over represented temperatures."""

    model: ThermalDeBroglieWavelengthModel
    temperatures: VectorQuantity
    wavelength_unit: PhysicalUnit

    def __post_init__(self) -> None:
        """Check each wavelength-evaluation request argument."""
        self._check_arg_model()
        self._check_arg_temperatures()
        self._check_arg_wavelength_unit()

    def _check_arg_model(self) -> None:
        """Require the thermal-wavelength model."""
        if not isinstance(self.model, ThermalDeBroglieWavelengthModel):
            raise TypeError("model must be ThermalDeBroglieWavelengthModel")

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
        converted = self.temperatures_in_kelvin
        if converted.magnitude.size == 0:
            raise ValueError("temperatures must be nonempty")
        if np.any(converted.magnitude <= 0.0):
            raise ValueError("temperatures must be above absolute zero")

    def _check_arg_wavelength_unit(self) -> None:
        """Require an explicit physical length output unit."""
        if not isinstance(self.wavelength_unit, PhysicalUnit):
            raise TypeError("wavelength_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.wavelength_unit,
            _METRE,
        ):
            raise ValueError("wavelength_unit must be compatible with length")

    @property
    def temperatures_in_kelvin(self) -> VectorQuantity:
        """Return represented absolute temperatures in kelvin."""
        return ModelSystemUnitConverter.convert_vector(
            quantity=self.temperatures,
            target=_KELVIN,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class ThermalDeBroglieWavelengthEvaluation(ResultsObject):
    """Correlate one request with represented thermal wavelengths."""

    request: ThermalDeBroglieWavelengthEvaluationRequest
    wavelengths: VectorQuantity

    def __post_init__(self) -> None:
        """Check each wavelength-evaluation response argument."""
        self._check_arg_request()
        self._check_arg_wavelengths()

    def _check_arg_request(self) -> None:
        """Require the complete wavelength request."""
        if not isinstance(
            self.request,
            ThermalDeBroglieWavelengthEvaluationRequest,
        ):
            raise TypeError(
                "request must be ThermalDeBroglieWavelengthEvaluationRequest"
            )

    def _check_arg_wavelengths(self) -> None:
        """Require one positive wavelength per represented temperature."""
        if not isinstance(self.wavelengths, VectorQuantity):
            raise TypeError("wavelengths must be VectorQuantity")
        if self.wavelengths.unit != self.request.wavelength_unit:
            raise ValueError("wavelengths must use the requested unit")
        if self.wavelengths.magnitude.shape != (
            self.request.temperatures.magnitude.shape
        ):
            raise ValueError("wavelengths must match temperatures")
        if np.any(self.wavelengths.magnitude <= 0.0):
            raise ValueError("wavelengths must be positive")


@dataclass(frozen=True, slots=True)
class ThermalDeBroglieWavelengthEvaluator:
    r"""Evaluate $\lambda_T=h/\sqrt{2\pi m k_B T}$."""

    def action(
        self,
        *,
        request: ThermalDeBroglieWavelengthEvaluationRequest,
    ) -> ThermalDeBroglieWavelengthEvaluation:
        """Return wavelengths in the explicitly requested length unit."""
        if not isinstance(
            request,
            ThermalDeBroglieWavelengthEvaluationRequest,
        ):
            raise TypeError(
                "request must be ThermalDeBroglieWavelengthEvaluationRequest"
            )
        temperature_kelvin = request.temperatures_in_kelvin.magnitude
        mass_kilogram = request.model.particle_mass_in_kilograms.magnitude
        wavelength_metres = VectorQuantity(
            magnitude=(
                SI.h
                / np.sqrt(2.0 * math.pi * mass_kilogram * SI.k_B * temperature_kelvin)
            ),
            unit=_METRE,
        )
        return ThermalDeBroglieWavelengthEvaluation(
            request=request,
            wavelengths=ModelSystemUnitConverter.convert_vector(
                quantity=wavelength_metres,
                target=request.wavelength_unit,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ThermalDegeneracyParameterEvaluationRequest(DataObject):
    r"""Request $n\lambda_T^3$ for one thermodynamic point."""

    model: ThermalDeBroglieWavelengthModel
    number_density: ScalarQuantity
    temperature: ScalarQuantity
    wavelength_unit: PhysicalUnit

    def __post_init__(self) -> None:
        """Check each degeneracy-parameter request argument."""
        self._check_arg_model()
        self._check_arg_number_density()
        self._check_arg_temperature()
        self._check_arg_wavelength_unit()

    def _check_arg_model(self) -> None:
        """Require the thermal-wavelength model."""
        if not isinstance(self.model, ThermalDeBroglieWavelengthModel):
            raise TypeError("model must be ThermalDeBroglieWavelengthModel")

    def _check_arg_number_density(self) -> None:
        """Require a nonnegative physical number density."""
        if not isinstance(self.number_density, ScalarQuantity):
            raise TypeError("number_density must be ScalarQuantity")
        if not isinstance(self.number_density.unit, PhysicalUnit):
            raise TypeError("number_density must use a physical unit")
        if not ModelSystemUnitConverter.compatible(
            self.number_density.unit,
            _NUMBER_DENSITY_SI,
        ):
            raise ValueError("number_density must use inverse-volume units")
        if self.number_density.magnitude < 0.0:
            raise ValueError("number_density must be nonnegative")

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

    def _check_arg_wavelength_unit(self) -> None:
        """Require an explicit physical length output unit."""
        if not isinstance(self.wavelength_unit, PhysicalUnit):
            raise TypeError("wavelength_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.wavelength_unit,
            _METRE,
        ):
            raise ValueError("wavelength_unit must be compatible with length")

    @property
    def temperature_in_kelvin(self) -> ScalarQuantity:
        """Return absolute temperature in kelvin."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.temperature,
            target=_KELVIN,
        )

    @property
    def number_density_in_si(self) -> ScalarQuantity:
        """Return number density in reciprocal cubic metres."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.number_density,
            target=_NUMBER_DENSITY_SI,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class ThermalDegeneracyParameterEvaluation(ResultsObject):
    """Correlate one request with wavelength and unitless degeneracy parameter."""

    request: ThermalDegeneracyParameterEvaluationRequest
    wavelength: ScalarQuantity
    degeneracy_parameter: ScalarQuantity

    def __post_init__(self) -> None:
        """Check each degeneracy-parameter response argument."""
        self._check_arg_request()
        self._check_arg_wavelength()
        self._check_arg_degeneracy_parameter()

    def _check_arg_request(self) -> None:
        """Require the complete degeneracy-parameter request."""
        if not isinstance(
            self.request,
            ThermalDegeneracyParameterEvaluationRequest,
        ):
            raise TypeError(
                "request must be ThermalDegeneracyParameterEvaluationRequest"
            )

    def _check_arg_wavelength(self) -> None:
        """Require a positive wavelength in the requested unit."""
        if not isinstance(self.wavelength, ScalarQuantity):
            raise TypeError("wavelength must be ScalarQuantity")
        if self.wavelength.unit != self.request.wavelength_unit:
            raise ValueError("wavelength must use the requested unit")
        if self.wavelength.magnitude <= 0.0:
            raise ValueError("wavelength must be positive")

    def _check_arg_degeneracy_parameter(self) -> None:
        """Require a nonnegative explicitly unitless parameter."""
        if not isinstance(self.degeneracy_parameter, ScalarQuantity):
            raise TypeError("degeneracy_parameter must be ScalarQuantity")
        if not isinstance(self.degeneracy_parameter.unit, Unitless):
            raise ValueError("degeneracy_parameter must be unitless")
        if self.degeneracy_parameter.magnitude < 0.0:
            raise ValueError("degeneracy_parameter must be nonnegative")


@dataclass(frozen=True, slots=True)
class ThermalDegeneracyParameterEvaluator:
    r"""Evaluate the unitless phase-space-density parameter $n\lambda_T^3$."""

    def action(
        self,
        *,
        request: ThermalDegeneracyParameterEvaluationRequest,
    ) -> ThermalDegeneracyParameterEvaluation:
        """Return wavelength and degeneracy parameter for one state point."""
        if not isinstance(
            request,
            ThermalDegeneracyParameterEvaluationRequest,
        ):
            raise TypeError(
                "request must be ThermalDegeneracyParameterEvaluationRequest"
            )
        wavelength_evaluation = request.model.evaluate_wavelengths(
            temperatures=VectorQuantity(
                magnitude=np.array(
                    [request.temperature.magnitude],
                    dtype=np.float64,
                ),
                unit=request.temperature.unit,
            ),
            wavelength_unit=request.wavelength_unit,
        )
        wavelength = ScalarQuantity(
            magnitude=float(wavelength_evaluation.wavelengths.magnitude[0]),
            unit=request.wavelength_unit,
        )
        wavelength_metres = ModelSystemUnitConverter.convert_scalar(
            quantity=wavelength,
            target=_METRE,
        )
        parameter = (
            request.number_density_in_si.magnitude * wavelength_metres.magnitude**3
        )
        return ThermalDegeneracyParameterEvaluation(
            request=request,
            wavelength=wavelength,
            degeneracy_parameter=ScalarQuantity(
                magnitude=float(parameter),
                unit=Unitless(),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ThermalDeBroglieWavelengthModel(
    ActionizedDataObject[
        ThermalDeBroglieWavelengthEvaluationRequest,
        ThermalDeBroglieWavelengthEvaluation,
    ]
):
    """Represent particle mass for canonical thermal-wavelength evaluation."""

    particle_mass: ScalarQuantity
    actionizer: ThermalDeBroglieWavelengthEvaluator = field(
        default_factory=ThermalDeBroglieWavelengthEvaluator,
        init=False,
        repr=False,
        compare=False,
    )
    degeneracy_actionizer: ThermalDegeneracyParameterEvaluator = field(
        default_factory=ThermalDegeneracyParameterEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the particle-mass argument."""
        self._check_arg_particle_mass()

    def _check_arg_particle_mass(self) -> None:
        """Require a positive physical mass."""
        if not isinstance(self.particle_mass, ScalarQuantity):
            raise TypeError("particle_mass must be ScalarQuantity")
        if not isinstance(self.particle_mass.unit, PhysicalUnit):
            raise TypeError("particle_mass must use a physical mass unit")
        if not ModelSystemUnitConverter.compatible(
            self.particle_mass.unit,
            _KILOGRAM,
        ):
            raise ValueError("particle_mass must use units compatible with mass")
        if self.particle_mass.magnitude <= 0.0:
            raise ValueError("particle_mass must be positive")

    @property
    def particle_mass_in_kilograms(self) -> ScalarQuantity:
        """Return particle mass in kilograms."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.particle_mass,
            target=_KILOGRAM,
        )

    def evaluate_wavelengths(
        self,
        *,
        temperatures: VectorQuantity,
        wavelength_unit: PhysicalUnit,
    ) -> ThermalDeBroglieWavelengthEvaluation:
        """Evaluate thermal wavelengths over explicit temperatures."""
        return self._respond(
            request=ThermalDeBroglieWavelengthEvaluationRequest(
                model=self,
                temperatures=temperatures,
                wavelength_unit=wavelength_unit,
            )
        )

    def evaluate_degeneracy_parameter(
        self,
        *,
        number_density: ScalarQuantity,
        temperature: ScalarQuantity,
        wavelength_unit: PhysicalUnit,
    ) -> ThermalDegeneracyParameterEvaluation:
        r"""Evaluate $n\lambda_T^3$ without applying a regime threshold."""
        return self.degeneracy_actionizer.action(
            request=ThermalDegeneracyParameterEvaluationRequest(
                model=self,
                number_density=number_density,
                temperature=temperature,
                wavelength_unit=wavelength_unit,
            )
        )
