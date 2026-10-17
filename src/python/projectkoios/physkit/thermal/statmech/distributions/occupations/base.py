r"""Shared unit-aware records for ideal occupation distributions."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np
import numpy.typing as npt

from projectkoios.physkit.constants import SI
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

_JOULE = PhysicalUnit(expression="joule")
_KELVIN = PhysicalUnit(expression="kelvin")


class IdealOccupationStatistic(Enum):
    """Identify one ideal occupation formula."""

    FERMI_DIRAC = "fermi_dirac"
    BOSE_EINSTEIN = "bose_einstein"
    MAXWELL_BOLTZMANN = "maxwell_boltzmann"


@dataclass(frozen=True, slots=True, kw_only=True)
class IdealOccupationEvaluationRequest(DataObject):
    """Request occupation factors on explicit physical energy coordinates."""

    statistic: IdealOccupationStatistic
    energies: VectorQuantity
    chemical_potential: ScalarQuantity
    temperature: ScalarQuantity

    def __post_init__(self) -> None:
        """Check each occupation-evaluation request argument."""
        self._check_arg_statistic()
        self._check_arg_energies()
        self._check_arg_chemical_potential()
        self._check_arg_temperature()

    def _check_arg_statistic(self) -> None:
        """Require one closed ideal-occupation statistic member."""
        if not isinstance(self.statistic, IdealOccupationStatistic):
            raise TypeError("statistic must be IdealOccupationStatistic")

    def _check_arg_energies(self) -> None:
        """Require nonempty physical energy coordinates."""
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if not isinstance(self.energies.unit, PhysicalUnit):
            raise TypeError("energies must use a physical energy unit")
        if not ModelSystemUnitConverter.compatible(self.energies.unit, _JOULE):
            raise ValueError("energies must use units compatible with energy")
        if self.energies.magnitude.size == 0:
            raise ValueError("energies must be nonempty")

    def _check_arg_chemical_potential(self) -> None:
        """Require a physical chemical-potential energy."""
        if not isinstance(self.chemical_potential, ScalarQuantity):
            raise TypeError("chemical_potential must be ScalarQuantity")
        if not isinstance(self.chemical_potential.unit, PhysicalUnit):
            raise TypeError("chemical_potential must use a physical energy unit")
        if not ModelSystemUnitConverter.compatible(
            self.chemical_potential.unit,
            _JOULE,
        ):
            raise ValueError("chemical_potential must use units compatible with energy")

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
        if self.statistic is IdealOccupationStatistic.BOSE_EINSTEIN and np.any(
            self.dimensionless_energy_offsets <= 0.0
        ):
            raise ValueError(
                "Bose-Einstein energies must exceed the chemical potential"
            )

    @property
    def energies_in_joule(self) -> VectorQuantity:
        """Return represented energies in joules."""
        return ModelSystemUnitConverter.convert_vector(
            quantity=self.energies,
            target=_JOULE,
        )

    @property
    def chemical_potential_in_joule(self) -> ScalarQuantity:
        """Return chemical potential in joules."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.chemical_potential,
            target=_JOULE,
        )

    @property
    def temperature_in_kelvin(self) -> ScalarQuantity:
        """Return absolute temperature in kelvin."""
        return ModelSystemUnitConverter.convert_scalar(
            quantity=self.temperature,
            target=_KELVIN,
        )

    @property
    def dimensionless_energy_offsets(self) -> npt.NDArray[np.float64]:
        r"""Return $(\varepsilon-\mu)/(k_BT)$ in binary64 arithmetic."""
        return np.asarray(
            (
                self.energies_in_joule.magnitude
                - self.chemical_potential_in_joule.magnitude
            )
            / (SI.k_B * self.temperature_in_kelvin.magnitude),
            dtype=np.float64,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class IdealOccupationEvaluation(ResultsObject):
    """Correlate one physical request with unitless occupation factors."""

    request: IdealOccupationEvaluationRequest
    occupation_factors: VectorQuantity

    def __post_init__(self) -> None:
        """Check each occupation-evaluation response argument."""
        self._check_arg_request()
        self._check_arg_occupation_factors()

    def _check_arg_request(self) -> None:
        """Require the complete occupation request."""
        if not isinstance(self.request, IdealOccupationEvaluationRequest):
            raise TypeError("request must be IdealOccupationEvaluationRequest")

    def _check_arg_occupation_factors(self) -> None:
        """Require one finite nonnegative unitless factor per energy."""
        if not isinstance(self.occupation_factors, VectorQuantity):
            raise TypeError("occupation_factors must be VectorQuantity")
        if not isinstance(self.occupation_factors.unit, Unitless):
            raise ValueError("occupation_factors must be unitless")
        if (
            self.occupation_factors.magnitude.shape
            != self.request.energies.magnitude.shape
        ):
            raise ValueError("occupation_factors must match energies")
        if np.any(self.occupation_factors.magnitude < 0.0):
            raise ValueError("occupation_factors must be nonnegative")
