"""Shared typed contracts for radial pair-potential evaluation."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    VectorQuantity,
)

_METRE = PhysicalUnit(expression="meter")
_JOULE = PhysicalUnit(expression="joule")
_NEWTON = PhysicalUnit(expression="newton")


class RadialPairPotential(DataObject, ABC):
    """Define the supported operation of one immutable radial pair model."""

    __slots__ = ()

    @abstractmethod
    def evaluate(
        self,
        *,
        distances: VectorQuantity,
        energy_unit: PhysicalUnit,
        force_unit: PhysicalUnit,
    ) -> RadialPairPotentialEvaluation:
        """Evaluate potential energy and signed radial force."""


@dataclass(frozen=True, slots=True, kw_only=True)
class RadialPairPotentialEvaluationRequest(DataObject):
    """Request sampled energies and forces from one radial pair model."""

    model: RadialPairPotential
    distances: VectorQuantity
    energy_unit: PhysicalUnit
    force_unit: PhysicalUnit

    def __post_init__(self) -> None:
        """Check every radial-potential evaluation request argument."""
        self._check_arg_model()
        self._check_arg_distances()
        self._check_arg_energy_unit()
        self._check_arg_force_unit()

    def _check_arg_model(self) -> None:
        """Require one concrete radial pair-potential model."""
        if not isinstance(self.model, RadialPairPotential):
            raise TypeError("model must be RadialPairPotential")

    def _check_arg_distances(self) -> None:
        """Require a nonempty vector of positive physical lengths."""
        if not isinstance(self.distances, VectorQuantity):
            raise TypeError("distances must be VectorQuantity")
        if not isinstance(self.distances.unit, PhysicalUnit):
            raise TypeError("distances must use a physical length unit")
        if not ModelSystemUnitConverter.compatible(
            self.distances.unit,
            _METRE,
        ):
            raise ValueError("distances must use units compatible with length")
        if self.distances.magnitude.size == 0:
            raise ValueError("distances must be nonempty")
        if np.any(self.distances.magnitude <= 0.0):
            raise ValueError("distances must be positive")

    def _check_arg_energy_unit(self) -> None:
        """Require a physical energy output unit."""
        if not isinstance(self.energy_unit, PhysicalUnit):
            raise TypeError("energy_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(self.energy_unit, _JOULE):
            raise ValueError("energy_unit must be compatible with energy")

    def _check_arg_force_unit(self) -> None:
        """Require a physical force output unit."""
        if not isinstance(self.force_unit, PhysicalUnit):
            raise TypeError("force_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(self.force_unit, _NEWTON):
            raise ValueError("force_unit must be compatible with force")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class RadialPairPotentialEvaluation(ResultsObject):
    """Retain a radial-potential request and immutable sampled outputs."""

    request: RadialPairPotentialEvaluationRequest
    potential_energies: VectorQuantity
    radial_forces: VectorQuantity

    def __post_init__(self) -> None:
        """Check every radial-potential evaluation response argument."""
        self._check_arg_request()
        self._check_arg_potential_energies()
        self._check_arg_radial_forces()

    def _check_arg_request(self) -> None:
        """Require the complete evaluation request."""
        if not isinstance(self.request, RadialPairPotentialEvaluationRequest):
            raise TypeError("request must be RadialPairPotentialEvaluationRequest")

    def _check_arg_potential_energies(self) -> None:
        """Require one energy in the requested unit per distance."""
        if not isinstance(self.potential_energies, VectorQuantity):
            raise TypeError("potential_energies must be VectorQuantity")
        if self.potential_energies.unit != self.request.energy_unit:
            raise ValueError("potential_energies must use the requested energy unit")
        if (
            self.potential_energies.magnitude.shape
            != self.request.distances.magnitude.shape
        ):
            raise ValueError("potential_energies must match distances")

    def _check_arg_radial_forces(self) -> None:
        """Require one force in the requested unit per distance."""
        if not isinstance(self.radial_forces, VectorQuantity):
            raise TypeError("radial_forces must be VectorQuantity")
        if self.radial_forces.unit != self.request.force_unit:
            raise ValueError("radial_forces must use the requested force unit")
        if self.radial_forces.magnitude.shape != self.request.distances.magnitude.shape:
            raise ValueError("radial_forces must match distances")
