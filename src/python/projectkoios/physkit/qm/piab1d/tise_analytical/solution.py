"""Closed-form energy solution for the one-dimensional particle in a box.

The energy formula and mode ordering are adapted from
``ksdft2effmass/analysis/model_systems/particle_in_box/model.py`` at commit
``7bd913151f7e61ed2bdba593df920be36573b502``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.numerics.typing.numpy.arrays import IntegerVector
from projectkoios.physkit.qm.eigenfunctions import (
    SampledEigenfunction1D,
    SampledEigenfunctions1D,
)
from projectkoios.physkit.units import (
    ComplexVectorQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

from ..base import Piab1D


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class Piab1DAnalyticalResults(ResultsObject):
    """Retain one model and its ordered closed-form energy levels."""

    model: Piab1D
    quantum_numbers: IntegerVector
    energies: VectorQuantity

    def __post_init__(self) -> None:
        """Check each analytical-spectrum response argument."""
        self._check_arg_model()
        self._check_arg_quantum_numbers()
        self._check_arg_energies()

    def _check_arg_model(self) -> None:
        """Require the PIAB1D model that owns the analytical spectrum."""
        if not isinstance(self.model, Piab1D):
            raise TypeError("model must be Piab1D")

    def _check_arg_quantum_numbers(self) -> None:
        """Require consecutive positive integer eigenstate labels."""
        if not isinstance(self.quantum_numbers, np.ndarray):
            raise TypeError("quantum_numbers must be a numpy.ndarray")
        values = np.asarray(self.quantum_numbers)
        if values.ndim != 1 or values.dtype.kind not in "iu":
            raise ValueError("quantum_numbers must be an integer vector")
        expected = np.arange(1, values.size + 1, dtype=np.int64)
        if values.size == 0 or not np.array_equal(values, expected):
            raise ValueError("quantum_numbers must be consecutive positive integers")
        immutable = np.frombuffer(
            values.astype("<i8").tobytes(),
            dtype="<i8",
        )
        object.__setattr__(self, "quantum_numbers", immutable)

    def _check_arg_energies(self) -> None:
        """Require one model-unit eigenvalue for every quantum number."""
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if self.energies.magnitude.shape != self.quantum_numbers.shape:
            raise ValueError("energies must match quantum_numbers")
        if self.energies.unit != self.model.unit_system.energy_unit:
            raise ValueError("energies must use the model energy unit")

    def eigenfunctions(
        self,
        *,
        coordinates: VectorQuantity,
    ) -> SampledEigenfunctions1D:
        """Evaluate all represented normalized stationary eigenvectors."""
        if not isinstance(coordinates, VectorQuantity):
            raise TypeError("coordinates must be VectorQuantity")
        selected = ModelSystemUnitConverter.convert_vector(
            coordinates,
            self.model.unit_system.length_unit,
        )
        length = self.model.length_quantity.magnitude
        if np.any(selected.magnitude < 0.0) or np.any(selected.magnitude > length):
            raise ValueError("coordinates must lie in the closed box interval")
        indices = self.quantum_numbers.astype(np.float64)
        values = np.sqrt(2.0 / length) * np.sin(
            np.pi * np.outer(selected.magnitude / length, indices)
        )
        length_unit = self.model.unit_system.length_unit
        if isinstance(length_unit, Unitless):
            amplitude_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(length_unit, PhysicalUnit):
                raise TypeError("physical model requires a physical length unit")
            amplitude_unit = PhysicalUnit(
                expression=f"({length_unit.expression}) ** -0.5"
            )
        eigenfunctions = tuple(
            SampledEigenfunction1D(
                eigenvalue=ScalarQuantity(
                    magnitude=float(self.energies.magnitude[index]),
                    unit=self.energies.unit,
                ),
                coordinates=coordinates,
                amplitudes=ComplexVectorQuantity(
                    magnitude=np.asarray(values[:, index], dtype=np.complex128),
                    unit=amplitude_unit,
                ),
            )
            for index in range(self.quantum_numbers.size)
        )
        return SampledEigenfunctions1D(eigenfunctions=eigenfunctions)


@dataclass(frozen=True, slots=True, kw_only=True)
class Piab1DAnalyticalSolution:
    """Evaluate the known continuum energy spectrum of one ``Piab1D``."""

    model: Piab1D

    def __post_init__(self) -> None:
        """Check the configured PIAB1D model argument."""
        self._check_arg_model()

    def _check_arg_model(self) -> None:
        """Require the PIAB1D model whose spectrum will be evaluated."""
        if not isinstance(self.model, Piab1D):
            raise TypeError("model must be Piab1D")

    def evaluate(
        self,
        *,
        count: int,
    ) -> Piab1DAnalyticalResults:
        """Return the first ``count`` closed-form energy levels."""
        if type(count) is not int:
            raise TypeError("count must be a built-in int")
        if count <= 0:
            raise ValueError("count must be positive")
        quantum_numbers = np.arange(1, count + 1, dtype=np.int64)
        indices = quantum_numbers.astype(np.float64)
        hbar = self.model.hbar_quantity
        mass = self.model.mass_quantity
        length = self.model.length_quantity
        energies = (
            hbar.magnitude**2
            * np.pi**2
            * indices**2
            / (2.0 * mass.magnitude * length.magnitude**2)
        )
        if isinstance(self.model.unit_system.energy_unit, Unitless):
            energy_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(hbar.unit, PhysicalUnit):
                raise TypeError("physical model requires physical action units")
            if not isinstance(mass.unit, PhysicalUnit):
                raise TypeError("physical model requires physical mass units")
            if not isinstance(length.unit, PhysicalUnit):
                raise TypeError("physical model requires physical length units")
            source_unit = PhysicalUnit(
                expression=f"({hbar.unit.expression}) ** 2 / "
                f"({mass.unit.expression}) / "
                f"({length.unit.expression}) ** 2"
            )
            energy_unit = self.model.unit_system.energy_unit
            energies *= ModelSystemUnitConverter.conversion_factor(
                source_unit,
                energy_unit,
            )
        return Piab1DAnalyticalResults(
            model=self.model,
            quantum_numbers=quantum_numbers,
            energies=VectorQuantity(
                magnitude=energies,
                unit=energy_unit,
            ),
        )
