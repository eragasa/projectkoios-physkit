r"""Unit-aware symmetric binary regular-solution mixing model."""

from __future__ import annotations

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
_MOLAR_ENERGY_SI = PhysicalUnit(expression="joule / mole")
_MOLAR_ENTROPY_SI = PhysicalUnit(expression="joule / mole / kelvin")


@dataclass(frozen=True, slots=True, kw_only=True)
class RegularSolutionMixingEvaluationRequest(DataObject):
    """Request symmetric regular-solution mixing quantities."""

    model: SymmetricRegularSolutionModel
    component_b_mole_fractions: VectorQuantity
    temperature: ScalarQuantity
    molar_energy_unit: PhysicalUnit

    def __post_init__(self) -> None:
        """Check each mixing-evaluation request argument."""
        self._check_arg_model()
        self._check_arg_component_b_mole_fractions()
        self._check_arg_temperature()
        self._check_arg_molar_energy_unit()

    def _check_arg_model(self) -> None:
        """Require the symmetric regular-solution model."""
        if not isinstance(self.model, SymmetricRegularSolutionModel):
            raise TypeError("model must be SymmetricRegularSolutionModel")

    def _check_arg_component_b_mole_fractions(self) -> None:
        """Require nonempty unitless binary mole fractions in $[0,1]$."""
        if not isinstance(self.component_b_mole_fractions, VectorQuantity):
            raise TypeError("component_b_mole_fractions must be VectorQuantity")
        if not isinstance(self.component_b_mole_fractions.unit, Unitless):
            raise ValueError("component_b_mole_fractions must be unitless")
        values = self.component_b_mole_fractions.magnitude
        if values.size == 0:
            raise ValueError("component_b_mole_fractions must be nonempty")
        if np.any(values < 0.0) or np.any(values > 1.0):
            raise ValueError(
                "component_b_mole_fractions must lie in the closed interval [0, 1]"
            )

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

    def _check_arg_molar_energy_unit(self) -> None:
        """Require an explicit molar-energy output unit."""
        if not isinstance(self.molar_energy_unit, PhysicalUnit):
            raise TypeError("molar_energy_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.molar_energy_unit,
            _MOLAR_ENERGY_SI,
        ):
            raise ValueError("molar_energy_unit must be compatible with molar energy")

    @property
    def temperature_in_kelvin(self) -> ScalarQuantity:
        """Return absolute temperature in kelvin."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.temperature,
            target=_KELVIN,
        )

    @property
    def molar_entropy_unit(self) -> PhysicalUnit:
        """Return the entropy unit associated with the requested energy unit."""
        return PhysicalUnit(
            expression=f"({self.molar_energy_unit.expression}) / kelvin"
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class RegularSolutionMixingEvaluation(ResultsObject):
    """Correlate a request with molar mixing enthalpy, entropy, and free energy."""

    request: RegularSolutionMixingEvaluationRequest
    molar_enthalpy_of_mixing: VectorQuantity
    molar_entropy_of_mixing: VectorQuantity
    molar_gibbs_free_energy_of_mixing: VectorQuantity

    def __post_init__(self) -> None:
        """Check each mixing-evaluation response argument."""
        self._check_arg_request()
        self._check_arg_molar_enthalpy_of_mixing()
        self._check_arg_molar_entropy_of_mixing()
        self._check_arg_molar_gibbs_free_energy_of_mixing()

    def _check_arg_request(self) -> None:
        """Require the complete mixing request."""
        if not isinstance(self.request, RegularSolutionMixingEvaluationRequest):
            raise TypeError("request must be RegularSolutionMixingEvaluationRequest")

    def _check_arg_molar_enthalpy_of_mixing(self) -> None:
        """Require correlated molar mixing enthalpy."""
        self._check_molar_vector(
            quantity=self.molar_enthalpy_of_mixing,
            expected_unit=self.request.molar_energy_unit,
            name="molar_enthalpy_of_mixing",
        )

    def _check_arg_molar_entropy_of_mixing(self) -> None:
        """Require correlated molar mixing entropy."""
        self._check_molar_vector(
            quantity=self.molar_entropy_of_mixing,
            expected_unit=self.request.molar_entropy_unit,
            name="molar_entropy_of_mixing",
        )

    def _check_arg_molar_gibbs_free_energy_of_mixing(self) -> None:
        """Require correlated molar Gibbs free energy of mixing."""
        self._check_molar_vector(
            quantity=self.molar_gibbs_free_energy_of_mixing,
            expected_unit=self.request.molar_energy_unit,
            name="molar_gibbs_free_energy_of_mixing",
        )

    def _check_molar_vector(
        self,
        *,
        quantity: VectorQuantity,
        expected_unit: PhysicalUnit,
        name: str,
    ) -> None:
        """Check shared mechanical invariants of one output vector."""
        if not isinstance(quantity, VectorQuantity):
            raise TypeError(f"{name} must be VectorQuantity")
        if quantity.unit != expected_unit:
            raise ValueError(f"{name} must use its requested unit")
        if quantity.magnitude.shape != (
            self.request.component_b_mole_fractions.magnitude.shape
        ):
            raise ValueError(f"{name} must match component_b_mole_fractions")


@dataclass(frozen=True, slots=True)
class RegularSolutionMixingEvaluator:
    """Evaluate symmetric regular-solution molar mixing quantities."""

    def action(
        self,
        *,
        request: RegularSolutionMixingEvaluationRequest,
    ) -> RegularSolutionMixingEvaluation:
        """Return enthalpy, entropy, and Gibbs free energy of mixing."""
        if not isinstance(request, RegularSolutionMixingEvaluationRequest):
            raise TypeError("request must be RegularSolutionMixingEvaluationRequest")
        mole_fraction = request.component_b_mole_fractions.magnitude
        interaction = request.model.interaction_parameter_in_si.magnitude
        temperature = request.temperature_in_kelvin.magnitude
        entropy_shape = np.zeros_like(mole_fraction)
        interior = (mole_fraction > 0.0) & (mole_fraction < 1.0)
        interior_fraction = mole_fraction[interior]
        entropy_shape[interior] = -(
            interior_fraction * np.log(interior_fraction)
            + (1.0 - interior_fraction) * np.log(1.0 - interior_fraction)
        )
        enthalpy_si = interaction * mole_fraction * (1.0 - mole_fraction)
        entropy_si = SI.R_g * entropy_shape
        gibbs_si = enthalpy_si - temperature * entropy_si
        return RegularSolutionMixingEvaluation(
            request=request,
            molar_enthalpy_of_mixing=self._convert_vector(
                values=enthalpy_si,
                source=_MOLAR_ENERGY_SI,
                target=request.molar_energy_unit,
            ),
            molar_entropy_of_mixing=self._convert_vector(
                values=entropy_si,
                source=_MOLAR_ENTROPY_SI,
                target=request.molar_entropy_unit,
            ),
            molar_gibbs_free_energy_of_mixing=self._convert_vector(
                values=gibbs_si,
                source=_MOLAR_ENERGY_SI,
                target=request.molar_energy_unit,
            ),
        )

    @staticmethod
    def _convert_vector(
        *,
        values: np.ndarray,
        source: PhysicalUnit,
        target: PhysicalUnit,
    ) -> VectorQuantity:
        """Convert one finite output vector to its requested unit."""
        return ModelSystemUnitConverter.convert_vector(
            quantity=VectorQuantity(
                magnitude=np.asarray(values, dtype=np.float64),
                unit=source,
            ),
            target=target,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class SymmetricRegularSolutionModel(
    ActionizedDataObject[
        RegularSolutionMixingEvaluationRequest,
        RegularSolutionMixingEvaluation,
    ]
):
    """Represent the composition-independent binary interaction parameter."""

    interaction_parameter: ScalarQuantity
    actionizer: RegularSolutionMixingEvaluator = field(
        default_factory=RegularSolutionMixingEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the interaction-parameter argument."""
        self._check_arg_interaction_parameter()

    def _check_arg_interaction_parameter(self) -> None:
        """Require a physical molar-energy interaction parameter."""
        if not isinstance(self.interaction_parameter, ScalarQuantity):
            raise TypeError("interaction_parameter must be ScalarQuantity")
        if not isinstance(self.interaction_parameter.unit, PhysicalUnit):
            raise TypeError("interaction_parameter must use a physical unit")
        if not ModelSystemUnitConverter.compatible(
            self.interaction_parameter.unit,
            _MOLAR_ENERGY_SI,
        ):
            raise ValueError(
                "interaction_parameter must use units compatible with molar energy"
            )

    @property
    def interaction_parameter_in_si(self) -> ScalarQuantity:
        """Return the interaction parameter in joules per mole."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.interaction_parameter,
            target=_MOLAR_ENERGY_SI,
        )

    def evaluate_mixing(
        self,
        *,
        component_b_mole_fractions: VectorQuantity,
        temperature: ScalarQuantity,
        molar_energy_unit: PhysicalUnit,
    ) -> RegularSolutionMixingEvaluation:
        """Evaluate molar mixing quantities over binary composition."""
        return self._respond(
            request=RegularSolutionMixingEvaluationRequest(
                model=self,
                component_b_mole_fractions=component_b_mole_fractions,
                temperature=temperature,
                molar_energy_unit=molar_energy_unit,
            )
        )
