r"""Unit-aware blackbody spectral energy density per wavelength."""

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
    VectorQuantity,
)

_KELVIN = PhysicalUnit(expression="kelvin")
_METRE = PhysicalUnit(expression="meter")
_SPECTRAL_ENERGY_DENSITY_SI = PhysicalUnit(expression="joule / meter ** 4")


@dataclass(frozen=True, slots=True, kw_only=True)
class BlackbodyWavelengthSpectrumEvaluationRequest(DataObject):
    """Request wavelength-domain blackbody spectra at one temperature."""

    model: BlackbodyWavelengthSpectrumModel
    wavelengths: VectorQuantity
    temperature: ScalarQuantity
    spectral_energy_density_unit: PhysicalUnit

    def __post_init__(self) -> None:
        """Check each spectrum-evaluation request argument."""
        self._check_arg_model()
        self._check_arg_wavelengths()
        self._check_arg_temperature()
        self._check_arg_spectral_energy_density_unit()

    def _check_arg_model(self) -> None:
        """Require the wavelength-spectrum model."""
        if not isinstance(self.model, BlackbodyWavelengthSpectrumModel):
            raise TypeError("model must be BlackbodyWavelengthSpectrumModel")

    def _check_arg_wavelengths(self) -> None:
        """Require nonempty positive physical wavelengths."""
        if not isinstance(self.wavelengths, VectorQuantity):
            raise TypeError("wavelengths must be VectorQuantity")
        if not isinstance(self.wavelengths.unit, PhysicalUnit):
            raise TypeError("wavelengths must use a physical length unit")
        if not ModelSystemUnitConverter.compatible(
            self.wavelengths.unit,
            _METRE,
        ):
            raise ValueError("wavelengths must use units compatible with length")
        if self.wavelengths.magnitude.size == 0:
            raise ValueError("wavelengths must be nonempty")
        if np.any(self.wavelengths.magnitude <= 0.0):
            raise ValueError("wavelengths must be positive")

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

    def _check_arg_spectral_energy_density_unit(self) -> None:
        """Require an explicit compatible spectral-density output unit."""
        if not isinstance(self.spectral_energy_density_unit, PhysicalUnit):
            raise TypeError("spectral_energy_density_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.spectral_energy_density_unit,
            _SPECTRAL_ENERGY_DENSITY_SI,
        ):
            raise ValueError(
                "spectral_energy_density_unit must be compatible with energy per "
                "volume per wavelength"
            )

    @property
    def wavelengths_in_metres(self) -> VectorQuantity:
        """Return represented wavelengths in metres."""
        return ModelSystemUnitConverter.convert_vector(
            quantity=self.wavelengths,
            target=_METRE,
        )

    @property
    def temperature_in_kelvin(self) -> ScalarQuantity:
        """Return absolute temperature in kelvin."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.temperature,
            target=_KELVIN,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class BlackbodyWavelengthSpectrumEvaluation(ResultsObject):
    """Correlate one request with exact and asymptotic spectral densities."""

    request: BlackbodyWavelengthSpectrumEvaluationRequest
    planck_spectral_energy_density: VectorQuantity
    wien_spectral_energy_density: VectorQuantity
    rayleigh_jeans_spectral_energy_density: VectorQuantity

    def __post_init__(self) -> None:
        """Check each spectrum-evaluation response argument."""
        self._check_arg_request()
        self._check_arg_planck_spectral_energy_density()
        self._check_arg_wien_spectral_energy_density()
        self._check_arg_rayleigh_jeans_spectral_energy_density()

    def _check_arg_request(self) -> None:
        """Require the complete spectrum request."""
        if not isinstance(
            self.request,
            BlackbodyWavelengthSpectrumEvaluationRequest,
        ):
            raise TypeError(
                "request must be BlackbodyWavelengthSpectrumEvaluationRequest"
            )

    def _check_arg_planck_spectral_energy_density(self) -> None:
        """Require a correlated nonnegative Planck spectrum."""
        self._check_spectral_energy_density(
            quantity=self.planck_spectral_energy_density,
            name="planck_spectral_energy_density",
        )

    def _check_arg_wien_spectral_energy_density(self) -> None:
        """Require a correlated nonnegative Wien spectrum."""
        self._check_spectral_energy_density(
            quantity=self.wien_spectral_energy_density,
            name="wien_spectral_energy_density",
        )

    def _check_arg_rayleigh_jeans_spectral_energy_density(self) -> None:
        """Require a correlated nonnegative Rayleigh--Jeans spectrum."""
        self._check_spectral_energy_density(
            quantity=self.rayleigh_jeans_spectral_energy_density,
            name="rayleigh_jeans_spectral_energy_density",
        )

    def _check_spectral_energy_density(
        self,
        *,
        quantity: VectorQuantity,
        name: str,
    ) -> None:
        """Check shared mechanical invariants of one returned spectrum."""
        if not isinstance(quantity, VectorQuantity):
            raise TypeError(f"{name} must be VectorQuantity")
        if quantity.unit != self.request.spectral_energy_density_unit:
            raise ValueError(f"{name} must use the requested unit")
        if quantity.magnitude.shape != self.request.wavelengths.magnitude.shape:
            raise ValueError(f"{name} must match wavelengths")
        if np.any(quantity.magnitude < 0.0):
            raise ValueError(f"{name} must be nonnegative")


@dataclass(frozen=True, slots=True)
class BlackbodyWavelengthSpectrumEvaluator:
    """Evaluate Planck density and its Wien and Rayleigh--Jeans limits."""

    def action(
        self,
        *,
        request: BlackbodyWavelengthSpectrumEvaluationRequest,
    ) -> BlackbodyWavelengthSpectrumEvaluation:
        """Return exact and asymptotic wavelength spectra."""
        if not isinstance(
            request,
            BlackbodyWavelengthSpectrumEvaluationRequest,
        ):
            raise TypeError(
                "request must be BlackbodyWavelengthSpectrumEvaluationRequest"
            )
        wavelength_metres = request.wavelengths_in_metres.magnitude
        temperature_kelvin = request.temperature_in_kelvin.magnitude
        exponent = SI.h * SI.c / (wavelength_metres * SI.k_B * temperature_kelvin)
        with np.errstate(over="ignore", under="ignore", divide="ignore"):
            exponential_factor = np.exp(-exponent)
            prefactor = 8.0 * math.pi * SI.h * SI.c / wavelength_metres**5
            planck_si = prefactor * exponential_factor / (-np.expm1(-exponent))
            wien_si = prefactor * exponential_factor
            rayleigh_jeans_si = (
                8.0 * math.pi * SI.k_B * temperature_kelvin / wavelength_metres**4
            )
        if not all(
            np.all(np.isfinite(values))
            for values in (planck_si, wien_si, rayleigh_jeans_si)
        ):
            raise ValueError(
                "wavelengths and temperature produce spectral densities outside "
                "the finite binary64 range"
            )
        return BlackbodyWavelengthSpectrumEvaluation(
            request=request,
            planck_spectral_energy_density=self._convert_density(
                values=planck_si,
                target=request.spectral_energy_density_unit,
            ),
            wien_spectral_energy_density=self._convert_density(
                values=wien_si,
                target=request.spectral_energy_density_unit,
            ),
            rayleigh_jeans_spectral_energy_density=self._convert_density(
                values=rayleigh_jeans_si,
                target=request.spectral_energy_density_unit,
            ),
        )

    @staticmethod
    def _convert_density(
        *,
        values: np.ndarray,
        target: PhysicalUnit,
    ) -> VectorQuantity:
        """Convert one SI spectral-density vector to the requested unit."""
        return ModelSystemUnitConverter.convert_vector(
            quantity=VectorQuantity(
                magnitude=np.asarray(values, dtype=np.float64),
                unit=_SPECTRAL_ENERGY_DENSITY_SI,
            ),
            target=target,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class BlackbodyWavelengthSpectrumModel(
    ActionizedDataObject[
        BlackbodyWavelengthSpectrumEvaluationRequest,
        BlackbodyWavelengthSpectrumEvaluation,
    ]
):
    """Expose unit-aware wavelength-domain blackbody spectrum evaluation."""

    actionizer: BlackbodyWavelengthSpectrumEvaluator = field(
        default_factory=BlackbodyWavelengthSpectrumEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def evaluate(
        self,
        *,
        wavelengths: VectorQuantity,
        temperature: ScalarQuantity,
        spectral_energy_density_unit: PhysicalUnit,
    ) -> BlackbodyWavelengthSpectrumEvaluation:
        """Evaluate exact and limiting spectra on physical wavelengths."""
        return self._respond(
            request=BlackbodyWavelengthSpectrumEvaluationRequest(
                model=self,
                wavelengths=wavelengths,
                temperature=temperature,
                spectral_energy_density_unit=spectral_energy_density_unit,
            )
        )
