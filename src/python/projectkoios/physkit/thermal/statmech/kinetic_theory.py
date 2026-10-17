r"""Unit-aware classical kinetic-theory states and distributions.

The implemented model represents the Maxwell--Boltzmann speed probability
density for a dilute classical gas in thermal equilibrium. Physical inputs and
results retain explicit units; evaluation converts compatible inputs through
the package's maintained unit converter.
"""

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
_KILOGRAM_PER_MOLE = PhysicalUnit(expression="kilogram / mole")
_METRE_PER_SECOND = PhysicalUnit(expression="meter / second")
_SECOND_PER_METRE = PhysicalUnit(expression="second / meter")


@dataclass(frozen=True, slots=True, kw_only=True)
class MaxwellBoltzmannSpeedDistributionEvaluationRequest(DataObject):
    """Request a speed distribution for one gas state and speed grid."""

    state: MaxwellBoltzmannGasState
    speeds: VectorQuantity

    def __post_init__(self) -> None:
        """Check each speed-distribution request argument."""
        self._check_arg_state()
        self._check_arg_speeds()

    def _check_arg_state(self) -> None:
        """Require the Maxwell--Boltzmann gas state."""
        if not isinstance(self.state, MaxwellBoltzmannGasState):
            raise TypeError("state must be MaxwellBoltzmannGasState")

    def _check_arg_speeds(self) -> None:
        """Require nonnegative physical speeds in compatible units."""
        if not isinstance(self.speeds, VectorQuantity):
            raise TypeError("speeds must be VectorQuantity")
        if not isinstance(self.speeds.unit, PhysicalUnit):
            raise TypeError("speeds must use a physical speed unit")
        if not ModelSystemUnitConverter.compatible(
            self.speeds.unit,
            _METRE_PER_SECOND,
        ):
            raise ValueError("speeds must use units compatible with speed")
        if np.any(self.speeds.magnitude < 0.0):
            raise ValueError("speeds must be nonnegative")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class MaxwellBoltzmannSpeedDistributionEvaluation(ResultsObject):
    """Correlate a speed-distribution request with probability density."""

    request: MaxwellBoltzmannSpeedDistributionEvaluationRequest
    probability_density: VectorQuantity

    def __post_init__(self) -> None:
        """Check each speed-distribution response argument."""
        self._check_arg_request()
        self._check_arg_probability_density()

    def _check_arg_request(self) -> None:
        """Require the complete speed-distribution request."""
        if not isinstance(
            self.request,
            MaxwellBoltzmannSpeedDistributionEvaluationRequest,
        ):
            raise TypeError(
                "request must be MaxwellBoltzmannSpeedDistributionEvaluationRequest"
            )

    def _check_arg_probability_density(self) -> None:
        """Require nonnegative density in reciprocal requested-speed units."""
        if not isinstance(self.probability_density, VectorQuantity):
            raise TypeError("probability_density must be VectorQuantity")
        if (
            self.probability_density.magnitude.shape
            != self.request.speeds.magnitude.shape
        ):
            raise ValueError("probability_density must match speeds")
        if np.any(self.probability_density.magnitude < 0.0):
            raise ValueError("probability_density must be nonnegative")
        expected_unit = PhysicalUnit(
            expression=f"({self.request.speeds.unit.expression}) ** -1"
        )
        if self.probability_density.unit != expected_unit:
            raise ValueError(
                "probability_density must use reciprocal requested-speed units"
            )


@dataclass(frozen=True, slots=True)
class MaxwellBoltzmannSpeedDistributionEvaluator:
    r"""Evaluate the Maxwell--Boltzmann speed probability density."""

    def action(
        self,
        *,
        request: MaxwellBoltzmannSpeedDistributionEvaluationRequest,
    ) -> MaxwellBoltzmannSpeedDistributionEvaluation:
        """Return probability density with respect to requested speed units."""
        if not isinstance(
            request,
            MaxwellBoltzmannSpeedDistributionEvaluationRequest,
        ):
            raise TypeError(
                "request must be MaxwellBoltzmannSpeedDistributionEvaluationRequest"
            )
        temperature_kelvin = request.state.temperature_in_kelvin.magnitude
        molar_mass_kg_per_mol = request.state.molar_mass_in_si.magnitude
        speeds_si = ModelSystemUnitConverter.convert_vector(
            quantity=request.speeds,
            target=_METRE_PER_SECOND,
        )
        mass_temperature_factor = molar_mass_kg_per_mol / (
            2.0 * SI.R_g * temperature_kelvin
        )
        prefactor = (4.0 / np.sqrt(np.pi)) * mass_temperature_factor**1.5
        squared_speeds = np.square(speeds_si.magnitude)
        density_si = VectorQuantity(
            magnitude=(
                prefactor
                * squared_speeds
                * np.exp(-mass_temperature_factor * squared_speeds)
            ),
            unit=_SECOND_PER_METRE,
        )
        density_unit = PhysicalUnit(
            expression=f"({request.speeds.unit.expression}) ** -1"
        )
        return MaxwellBoltzmannSpeedDistributionEvaluation(
            request=request,
            probability_density=ModelSystemUnitConverter.convert_vector(
                quantity=density_si,
                target=density_unit,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class MaxwellBoltzmannGasState(
    ActionizedDataObject[
        MaxwellBoltzmannSpeedDistributionEvaluationRequest,
        MaxwellBoltzmannSpeedDistributionEvaluation,
    ]
):
    """Represent temperature and molar mass for a classical ideal gas."""

    temperature: ScalarQuantity
    molar_mass: ScalarQuantity
    actionizer: MaxwellBoltzmannSpeedDistributionEvaluator = field(
        default_factory=MaxwellBoltzmannSpeedDistributionEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check each Maxwell--Boltzmann gas-state argument."""
        self._check_arg_temperature()
        self._check_arg_molar_mass()

    def _check_arg_temperature(self) -> None:
        """Require a positive physical thermodynamic temperature."""
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

    def _check_arg_molar_mass(self) -> None:
        """Require a positive physical molar mass."""
        if not isinstance(self.molar_mass, ScalarQuantity):
            raise TypeError("molar_mass must be ScalarQuantity")
        if not isinstance(self.molar_mass.unit, PhysicalUnit):
            raise TypeError("molar_mass must use a physical molar-mass unit")
        if not ModelSystemUnitConverter.compatible(
            self.molar_mass.unit,
            _KILOGRAM_PER_MOLE,
        ):
            raise ValueError("molar_mass must use units compatible with molar mass")
        if self.molar_mass_in_si.magnitude <= 0.0:
            raise ValueError("molar_mass must be positive")

    @property
    def temperature_in_kelvin(self) -> ScalarQuantity:
        """Return absolute temperature in kelvin, including affine conversion."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.temperature,
            target=_KELVIN,
        )

    @property
    def molar_mass_in_si(self) -> ScalarQuantity:
        """Return molar mass in kilograms per mole."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.molar_mass,
            target=_KILOGRAM_PER_MOLE,
        )

    def _characteristic_speed(self, *, coefficient: float) -> ScalarQuantity:
        """Return one SI characteristic speed for a dimensionless coefficient."""
        temperature_kelvin = self.temperature_in_kelvin.magnitude
        molar_mass_kg_per_mol = self.molar_mass_in_si.magnitude
        return ScalarQuantity(
            magnitude=math.sqrt(
                coefficient * SI.R_g * temperature_kelvin / molar_mass_kg_per_mol
            ),
            unit=_METRE_PER_SECOND,
        )

    @property
    def most_probable_speed(self) -> ScalarQuantity:
        """Return the analytic most probable speed."""
        return self._characteristic_speed(coefficient=2.0)

    @property
    def mean_speed(self) -> ScalarQuantity:
        """Return the analytic mean speed."""
        return self._characteristic_speed(coefficient=8.0 / math.pi)

    @property
    def root_mean_square_speed(self) -> ScalarQuantity:
        """Return the analytic root-mean-square speed."""
        return self._characteristic_speed(coefficient=3.0)

    def evaluate_speed_distribution(
        self,
        *,
        speeds: VectorQuantity,
    ) -> MaxwellBoltzmannSpeedDistributionEvaluation:
        """Evaluate speed probability density on explicit physical speeds."""
        return self._respond(
            request=MaxwellBoltzmannSpeedDistributionEvaluationRequest(
                state=self,
                speeds=speeds,
            )
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class HardSphereIdealGasState(DataObject):
    """Represent an equilibrium identical-particle hard-sphere gas state."""

    pressure: ScalarQuantity
    temperature: ScalarQuantity
    collision_diameter: ScalarQuantity

    def __post_init__(self) -> None:
        """Check each hard-sphere ideal-gas state argument."""
        self._check_arg_pressure()
        self._check_arg_temperature()
        self._check_arg_collision_diameter()

    def _check_arg_pressure(self) -> None:
        """Require a positive physical pressure."""
        if not isinstance(self.pressure, ScalarQuantity):
            raise TypeError("pressure must be ScalarQuantity")
        if not isinstance(self.pressure.unit, PhysicalUnit):
            raise TypeError("pressure must use a physical pressure unit")
        pressure_unit = PhysicalUnit(expression="pascal")
        if not ModelSystemUnitConverter.compatible(
            self.pressure.unit,
            pressure_unit,
        ):
            raise ValueError("pressure must use units compatible with pressure")
        if self.pressure_in_pascal.magnitude <= 0.0:
            raise ValueError("pressure must be positive")

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

    def _check_arg_collision_diameter(self) -> None:
        """Require a positive physical length for the collision diameter."""
        if not isinstance(self.collision_diameter, ScalarQuantity):
            raise TypeError("collision_diameter must be ScalarQuantity")
        if not isinstance(self.collision_diameter.unit, PhysicalUnit):
            raise TypeError("collision_diameter must use a physical length unit")
        length_unit = PhysicalUnit(expression="meter")
        if not ModelSystemUnitConverter.compatible(
            self.collision_diameter.unit,
            length_unit,
        ):
            raise ValueError("collision_diameter must use units compatible with length")
        if self.collision_diameter_in_metres.magnitude <= 0.0:
            raise ValueError("collision_diameter must be positive")

    @property
    def pressure_in_pascal(self) -> ScalarQuantity:
        """Return pressure in pascals."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.pressure,
            target=PhysicalUnit(expression="pascal"),
        )

    @property
    def temperature_in_kelvin(self) -> ScalarQuantity:
        """Return absolute temperature in kelvin, including affine conversion."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.temperature,
            target=_KELVIN,
        )

    @property
    def collision_diameter_in_metres(self) -> ScalarQuantity:
        """Return collision diameter in metres."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.collision_diameter,
            target=PhysicalUnit(expression="meter"),
        )

    @property
    def collision_cross_section(self) -> ScalarQuantity:
        r"""Return the hard-sphere cross section $\sigma=\pi d^2$."""
        diameter_metres = self.collision_diameter_in_metres.magnitude
        return ScalarQuantity(
            magnitude=math.pi * diameter_metres**2,
            unit=PhysicalUnit(expression="meter ** 2"),
        )

    @property
    def number_density(self) -> ScalarQuantity:
        r"""Return ideal-gas number density $n=P/(k_B T)$."""
        return ScalarQuantity(
            magnitude=(
                self.pressure_in_pascal.magnitude
                / (SI.k_B * self.temperature_in_kelvin.magnitude)
            ),
            unit=PhysicalUnit(expression="meter ** -3"),
        )

    @property
    def mean_free_path(self) -> ScalarQuantity:
        r"""Return $\lambda=(\sqrt{2}\,n\sigma)^{-1}$ for identical particles."""
        return ScalarQuantity(
            magnitude=(
                1.0
                / (
                    math.sqrt(2.0)
                    * self.number_density.magnitude
                    * self.collision_cross_section.magnitude
                )
            ),
            unit=PhysicalUnit(expression="meter"),
        )
