r"""Signed planar-interface Hertz--Knudsen fluxes for a 3D ideal gas."""

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
_PASCAL = PhysicalUnit(expression="pascal")
_NUMBER_FLUX_SI = PhysicalUnit(expression="meter ** -2 / second")
_MASS_FLUX_SI = PhysicalUnit(expression="kilogram / meter ** 2 / second")


@dataclass(frozen=True, slots=True, kw_only=True)
class PlanarHertzKnudsen3DFluxEvaluationRequest(DataObject):
    """Request signed fluxes on correlated temperature and pressure vectors."""

    model: PlanarHertzKnudsen3DNetFluxModel
    temperatures: VectorQuantity
    equilibrium_pressures: VectorQuantity
    ambient_pressures: VectorQuantity
    number_flux_unit: PhysicalUnit
    mass_flux_unit: PhysicalUnit

    def __post_init__(self) -> None:
        """Check each Hertz--Knudsen request argument."""
        self._check_arg_model()
        self._check_arg_temperatures()
        self._check_arg_equilibrium_pressures()
        self._check_arg_ambient_pressures()
        self._check_arg_number_flux_unit()
        self._check_arg_mass_flux_unit()

    def _check_arg_model(self) -> None:
        """Require the signed net-flux model."""
        if not isinstance(self.model, PlanarHertzKnudsen3DNetFluxModel):
            raise TypeError("model must be PlanarHertzKnudsen3DNetFluxModel")

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
        if np.any(self.temperatures_in_kelvin.magnitude <= 0.0):
            raise ValueError("temperatures must be above absolute zero")

    def _check_arg_equilibrium_pressures(self) -> None:
        """Require nonnegative physical equilibrium pressures."""
        self._check_pressure_vector(
            quantity=self.equilibrium_pressures,
            name="equilibrium_pressures",
        )

    def _check_arg_ambient_pressures(self) -> None:
        """Require nonnegative physical ambient pressures."""
        self._check_pressure_vector(
            quantity=self.ambient_pressures,
            name="ambient_pressures",
        )

    def _check_arg_number_flux_unit(self) -> None:
        """Require an explicit compatible number-flux unit."""
        if not isinstance(self.number_flux_unit, PhysicalUnit):
            raise TypeError("number_flux_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.number_flux_unit,
            _NUMBER_FLUX_SI,
        ):
            raise ValueError("number_flux_unit must be compatible with number flux")

    def _check_arg_mass_flux_unit(self) -> None:
        """Require an explicit compatible mass-flux unit."""
        if not isinstance(self.mass_flux_unit, PhysicalUnit):
            raise TypeError("mass_flux_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.mass_flux_unit,
            _MASS_FLUX_SI,
        ):
            raise ValueError("mass_flux_unit must be compatible with mass flux")

    def _check_pressure_vector(
        self,
        *,
        quantity: VectorQuantity,
        name: str,
    ) -> None:
        """Check one pressure vector against shared request coordinates."""
        if not isinstance(quantity, VectorQuantity):
            raise TypeError(f"{name} must be VectorQuantity")
        if not isinstance(quantity.unit, PhysicalUnit):
            raise TypeError(f"{name} must use a physical pressure unit")
        if not ModelSystemUnitConverter.compatible(quantity.unit, _PASCAL):
            raise ValueError(f"{name} must use units compatible with pressure")
        if quantity.magnitude.shape != self.temperatures.magnitude.shape:
            raise ValueError(f"{name} must match temperatures")
        if np.any(quantity.magnitude < 0.0):
            raise ValueError(f"{name} must be nonnegative")

    @property
    def temperatures_in_kelvin(self) -> VectorQuantity:
        """Return represented temperatures in kelvin."""
        return ModelSystemUnitConverter.convert_vector(
            quantity=self.temperatures,
            target=_KELVIN,
        )

    @property
    def equilibrium_pressures_in_pascal(self) -> VectorQuantity:
        """Return equilibrium pressures in pascals."""
        return ModelSystemUnitConverter.convert_vector(
            quantity=self.equilibrium_pressures,
            target=_PASCAL,
        )

    @property
    def ambient_pressures_in_pascal(self) -> VectorQuantity:
        """Return ambient pressures in pascals."""
        return ModelSystemUnitConverter.convert_vector(
            quantity=self.ambient_pressures,
            target=_PASCAL,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PlanarHertzKnudsen3DFluxEvaluation(ResultsObject):
    """Correlate one request with evaporation, condensation, and signed fluxes."""

    request: PlanarHertzKnudsen3DFluxEvaluationRequest
    evaporation_number_flux: VectorQuantity
    condensation_number_flux: VectorQuantity
    net_number_flux: VectorQuantity
    net_mass_flux: VectorQuantity

    def __post_init__(self) -> None:
        """Check each Hertz--Knudsen response argument."""
        self._check_arg_request()
        self._check_arg_evaporation_number_flux()
        self._check_arg_condensation_number_flux()
        self._check_arg_net_number_flux()
        self._check_arg_net_mass_flux()

    def _check_arg_request(self) -> None:
        """Require the complete flux request."""
        if not isinstance(self.request, PlanarHertzKnudsen3DFluxEvaluationRequest):
            raise TypeError("request must be PlanarHertzKnudsen3DFluxEvaluationRequest")

    def _check_arg_evaporation_number_flux(self) -> None:
        """Require correlated nonnegative evaporation number flux."""
        self._check_flux_vector(
            quantity=self.evaporation_number_flux,
            expected_unit=self.request.number_flux_unit,
            name="evaporation_number_flux",
            require_nonnegative=True,
        )

    def _check_arg_condensation_number_flux(self) -> None:
        """Require correlated nonnegative condensation number flux."""
        self._check_flux_vector(
            quantity=self.condensation_number_flux,
            expected_unit=self.request.number_flux_unit,
            name="condensation_number_flux",
            require_nonnegative=True,
        )

    def _check_arg_net_number_flux(self) -> None:
        """Require correlated signed net number flux."""
        self._check_flux_vector(
            quantity=self.net_number_flux,
            expected_unit=self.request.number_flux_unit,
            name="net_number_flux",
            require_nonnegative=False,
        )

    def _check_arg_net_mass_flux(self) -> None:
        """Require correlated signed net mass flux."""
        self._check_flux_vector(
            quantity=self.net_mass_flux,
            expected_unit=self.request.mass_flux_unit,
            name="net_mass_flux",
            require_nonnegative=False,
        )

    def _check_flux_vector(
        self,
        *,
        quantity: VectorQuantity,
        expected_unit: PhysicalUnit,
        name: str,
        require_nonnegative: bool,
    ) -> None:
        """Check shared mechanical invariants of one returned flux."""
        if not isinstance(quantity, VectorQuantity):
            raise TypeError(f"{name} must be VectorQuantity")
        if quantity.unit != expected_unit:
            raise ValueError(f"{name} must use its requested unit")
        if quantity.magnitude.shape != self.request.temperatures.magnitude.shape:
            raise ValueError(f"{name} must match temperatures")
        if require_nonnegative and np.any(quantity.magnitude < 0.0):
            raise ValueError(f"{name} must be nonnegative")


@dataclass(frozen=True, slots=True)
class PlanarHertzKnudsen3DFluxEvaluator:
    """Evaluate separate 3D-gas contributions across a planar interface."""

    def action(
        self,
        *,
        request: PlanarHertzKnudsen3DFluxEvaluationRequest,
    ) -> PlanarHertzKnudsen3DFluxEvaluation:
        """Return component and signed net fluxes without clamping."""
        if not isinstance(request, PlanarHertzKnudsen3DFluxEvaluationRequest):
            raise TypeError("request must be PlanarHertzKnudsen3DFluxEvaluationRequest")
        temperature = request.temperatures_in_kelvin.magnitude
        particle_mass = request.model.particle_mass_in_kilograms.magnitude
        denominator = np.sqrt(2.0 * math.pi * particle_mass * SI.k_B * temperature)
        evaporation_si = (
            request.model.evaporation_coefficient.magnitude
            * request.equilibrium_pressures_in_pascal.magnitude
            / denominator
        )
        condensation_si = (
            request.model.condensation_coefficient.magnitude
            * request.ambient_pressures_in_pascal.magnitude
            / denominator
        )
        net_number_si = evaporation_si - condensation_si
        net_mass_si = particle_mass * net_number_si
        return PlanarHertzKnudsen3DFluxEvaluation(
            request=request,
            evaporation_number_flux=self._convert_flux(
                values=evaporation_si,
                source=_NUMBER_FLUX_SI,
                target=request.number_flux_unit,
            ),
            condensation_number_flux=self._convert_flux(
                values=condensation_si,
                source=_NUMBER_FLUX_SI,
                target=request.number_flux_unit,
            ),
            net_number_flux=self._convert_flux(
                values=net_number_si,
                source=_NUMBER_FLUX_SI,
                target=request.number_flux_unit,
            ),
            net_mass_flux=self._convert_flux(
                values=net_mass_si,
                source=_MASS_FLUX_SI,
                target=request.mass_flux_unit,
            ),
        )

    @staticmethod
    def _convert_flux(
        *,
        values: np.ndarray,
        source: PhysicalUnit,
        target: PhysicalUnit,
    ) -> VectorQuantity:
        """Convert one finite flux vector to its requested unit."""
        return ModelSystemUnitConverter.convert_vector(
            quantity=VectorQuantity(
                magnitude=np.asarray(values, dtype=np.float64),
                unit=source,
            ),
            target=target,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PlanarHertzKnudsen3DNetFluxModel(
    ActionizedDataObject[
        PlanarHertzKnudsen3DFluxEvaluationRequest,
        PlanarHertzKnudsen3DFluxEvaluation,
    ]
):
    """Represent a 3D ideal gas at a planar 2D phase boundary."""

    particle_mass: ScalarQuantity
    evaporation_coefficient: ScalarQuantity
    condensation_coefficient: ScalarQuantity
    actionizer: PlanarHertzKnudsen3DFluxEvaluator = field(
        default_factory=PlanarHertzKnudsen3DFluxEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check each signed net-flux model argument."""
        self._check_arg_particle_mass()
        self._check_arg_evaporation_coefficient()
        self._check_arg_condensation_coefficient()

    def _check_arg_particle_mass(self) -> None:
        """Require a positive physical particle mass."""
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

    def _check_arg_evaporation_coefficient(self) -> None:
        """Require an explicitly unitless coefficient in $[0,1]$."""
        self._check_coefficient(
            coefficient=self.evaporation_coefficient,
            name="evaporation_coefficient",
        )

    def _check_arg_condensation_coefficient(self) -> None:
        """Require an explicitly unitless coefficient in $[0,1]$."""
        self._check_coefficient(
            coefficient=self.condensation_coefficient,
            name="condensation_coefficient",
        )

    def _check_coefficient(
        self,
        *,
        coefficient: ScalarQuantity,
        name: str,
    ) -> None:
        """Check one interfacial transfer coefficient."""
        if not isinstance(coefficient, ScalarQuantity):
            raise TypeError(f"{name} must be ScalarQuantity")
        if not isinstance(coefficient.unit, Unitless):
            raise ValueError(f"{name} must be unitless")
        if coefficient.magnitude < 0.0 or coefficient.magnitude > 1.0:
            raise ValueError(f"{name} must lie in the closed interval [0, 1]")

    @property
    def particle_mass_in_kilograms(self) -> ScalarQuantity:
        """Return particle mass in kilograms."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.particle_mass,
            target=_KILOGRAM,
        )

    def evaluate(
        self,
        *,
        temperatures: VectorQuantity,
        equilibrium_pressures: VectorQuantity,
        ambient_pressures: VectorQuantity,
        number_flux_unit: PhysicalUnit,
        mass_flux_unit: PhysicalUnit,
    ) -> PlanarHertzKnudsen3DFluxEvaluation:
        """Evaluate signed interfacial number and mass fluxes."""
        return self._respond(
            request=PlanarHertzKnudsen3DFluxEvaluationRequest(
                model=self,
                temperatures=temperatures,
                equilibrium_pressures=equilibrium_pressures,
                ambient_pressures=ambient_pressures,
                number_flux_unit=number_flux_unit,
                mass_flux_unit=mass_flux_unit,
            )
        )
