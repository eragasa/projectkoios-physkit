r"""Normalized exponential Boltzmann energy distribution.

The represented continuous density assumes a constant density of states on
nonnegative energy. It is not the three-dimensional translational kinetic-energy
distribution, whose density-of-states factor is proportional to ``sqrt(E)``.
"""

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
    VectorQuantity,
)

_JOULE = PhysicalUnit(expression="joule")


@dataclass(frozen=True, slots=True, kw_only=True)
class BoltzmannExponentialEnergyEvaluationRequest(DataObject):
    """Request exponential Boltzmann density on represented energies."""

    distribution: BoltzmannExponentialEnergyDistribution
    energies: VectorQuantity

    def __post_init__(self) -> None:
        """Check each energy-density request argument."""
        self._check_arg_distribution()
        self._check_arg_energies()

    def _check_arg_distribution(self) -> None:
        """Require the normalized exponential distribution model."""
        if not isinstance(
            self.distribution,
            BoltzmannExponentialEnergyDistribution,
        ):
            raise TypeError(
                "distribution must be BoltzmannExponentialEnergyDistribution"
            )

    def _check_arg_energies(self) -> None:
        """Require nonempty nonnegative physical energies."""
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if not isinstance(self.energies.unit, PhysicalUnit):
            raise TypeError("energies must use a physical energy unit")
        if not ModelSystemUnitConverter.compatible(
            self.energies.unit,
            _JOULE,
        ):
            raise ValueError("energies must use units compatible with energy")
        if self.energies.magnitude.size == 0:
            raise ValueError("energies must be nonempty")
        if np.any(self.energies.magnitude < 0.0):
            raise ValueError("energies must be nonnegative")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class BoltzmannExponentialEnergyEvaluation(ResultsObject):
    """Correlate an energy request with normalized probability density."""

    request: BoltzmannExponentialEnergyEvaluationRequest
    probability_density: VectorQuantity

    def __post_init__(self) -> None:
        """Check each energy-density response argument."""
        self._check_arg_request()
        self._check_arg_probability_density()

    def _check_arg_request(self) -> None:
        """Require the complete energy-density request."""
        if not isinstance(
            self.request,
            BoltzmannExponentialEnergyEvaluationRequest,
        ):
            raise TypeError(
                "request must be BoltzmannExponentialEnergyEvaluationRequest"
            )

    def _check_arg_probability_density(self) -> None:
        """Require nonnegative density in reciprocal represented-energy units."""
        if not isinstance(self.probability_density, VectorQuantity):
            raise TypeError("probability_density must be VectorQuantity")
        if self.probability_density.magnitude.shape != (
            self.request.energies.magnitude.shape
        ):
            raise ValueError("probability_density must match energies")
        expected_unit = PhysicalUnit(
            expression=f"({self.request.energies.unit.expression}) ** -1"
        )
        if self.probability_density.unit != expected_unit:
            raise ValueError(
                "probability_density must use reciprocal represented-energy units"
            )
        if np.any(self.probability_density.magnitude < 0.0):
            raise ValueError("probability_density must be nonnegative")


@dataclass(frozen=True, slots=True)
class BoltzmannExponentialEnergyEvaluator:
    r"""Evaluate $f(E)=\exp(-E/\theta)/\theta$ for $E\geq0$."""

    def action(
        self,
        *,
        request: BoltzmannExponentialEnergyEvaluationRequest,
    ) -> BoltzmannExponentialEnergyEvaluation:
        """Return normalized density in reciprocal requested-energy units."""
        if not isinstance(
            request,
            BoltzmannExponentialEnergyEvaluationRequest,
        ):
            raise TypeError(
                "request must be BoltzmannExponentialEnergyEvaluationRequest"
            )
        thermal_energy = ModelSystemUnitConverter.convert_scalar(
            quantity=request.distribution.thermal_energy,
            target=request.energies.unit,
        )
        density_unit = PhysicalUnit(
            expression=f"({request.energies.unit.expression}) ** -1"
        )
        return BoltzmannExponentialEnergyEvaluation(
            request=request,
            probability_density=VectorQuantity(
                magnitude=(
                    np.exp(-request.energies.magnitude / thermal_energy.magnitude)
                    / thermal_energy.magnitude
                ),
                unit=density_unit,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class BoltzmannExponentialEnergyDistribution(
    ActionizedDataObject[
        BoltzmannExponentialEnergyEvaluationRequest,
        BoltzmannExponentialEnergyEvaluation,
    ]
):
    """Represent a normalized exponential density with thermal scale."""

    thermal_energy: ScalarQuantity
    actionizer: BoltzmannExponentialEnergyEvaluator = field(
        default_factory=BoltzmannExponentialEnergyEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the thermal-energy argument."""
        self._check_arg_thermal_energy()

    def _check_arg_thermal_energy(self) -> None:
        """Require a positive scalar physical energy."""
        if not isinstance(self.thermal_energy, ScalarQuantity):
            raise TypeError("thermal_energy must be ScalarQuantity")
        if not isinstance(self.thermal_energy.unit, PhysicalUnit):
            raise TypeError("thermal_energy must use a physical energy unit")
        if not ModelSystemUnitConverter.compatible(
            self.thermal_energy.unit,
            _JOULE,
        ):
            raise ValueError("thermal_energy must use units compatible with energy")
        if self.thermal_energy.magnitude <= 0.0:
            raise ValueError("thermal_energy must be positive")

    def evaluate(
        self,
        *,
        energies: VectorQuantity,
    ) -> BoltzmannExponentialEnergyEvaluation:
        """Evaluate density on explicit nonnegative physical energies."""
        return self._respond(
            request=BoltzmannExponentialEnergyEvaluationRequest(
                distribution=self,
                energies=energies,
            )
        )
