"""Closed-form spectrum of the three-dimensional rectangular particle in a box."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    Unitless,
    VectorQuantity,
)

from ...base import Piab3D

type QuantumNumberTriples = npt.NDArray[np.int64]


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class Piab3DAnalyticalSolution(ResultsObject):
    """Retain one model and its ordered closed-form stationary-state energies."""

    model: Piab3D
    quantum_numbers: QuantumNumberTriples
    energies: VectorQuantity

    def __post_init__(self) -> None:
        """Check each analytical-spectrum response argument."""
        self._check_arg_model()
        self._check_arg_quantum_numbers()
        self._check_arg_energies()

    def _check_arg_model(self) -> None:
        """Require the cuboid PIAB3D model that owns the spectrum."""
        if not isinstance(self.model, Piab3D):
            raise TypeError("model must be Piab3D")

    def _check_arg_quantum_numbers(self) -> None:
        """Require unique positive triples in deterministic within-tie order."""
        if not isinstance(self.quantum_numbers, np.ndarray):
            raise TypeError("quantum_numbers must be a numpy.ndarray")
        values = np.asarray(self.quantum_numbers)
        if values.ndim != 2 or values.shape[1:] != (3,):
            raise ValueError("quantum_numbers must have shape (state_count, 3)")
        if values.dtype.kind not in "iu":
            raise ValueError("quantum_numbers must contain integers")
        if values.shape[0] == 0:
            raise ValueError("quantum_numbers must be nonempty")
        if np.any(values <= 0):
            raise ValueError("quantum_numbers must be positive")
        if np.unique(values, axis=0).shape[0] != values.shape[0]:
            raise ValueError("quantum_numbers must be unique")
        immutable = np.frombuffer(values.astype("<i8").tobytes(), dtype="<i8")
        object.__setattr__(
            self,
            "quantum_numbers",
            immutable.reshape(values.shape),
        )

    def _check_arg_energies(self) -> None:
        """Require model-unit eigenvalues ordered with their eigenstate triples."""
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        values = self.quantum_numbers
        if self.energies.magnitude.shape != (values.shape[0],):
            raise ValueError("energies must match quantum_numbers")
        if self.energies.unit != self.model.unit_system.energy_unit:
            raise ValueError("energies must use the model energy unit")
        if np.any(self.energies.magnitude[1:] < self.energies.magnitude[:-1]):
            raise ValueError("energies must be nondecreasing")
        for index in range(1, values.shape[0]):
            if self.energies.magnitude[index] == self.energies.magnitude[
                index - 1
            ] and tuple(values[index]) < tuple(values[index - 1]):
                raise ValueError(
                    "equal-energy quantum_numbers must be lexicographically ordered"
                )


@dataclass(frozen=True, slots=True)
class Piab3DAnalyticalEvaluator:
    """Evaluate a finite rectangular-cuboid inventory of continuum energies."""

    def execute(
        self,
        *,
        model: Piab3D,
        maximum_quantum_number_x: int,
        maximum_quantum_number_y: int,
        maximum_quantum_number_z: int,
    ) -> Piab3DAnalyticalSolution:
        r"""Return ordered energies for selected positive quantum-number triples.

        The represented energies are

        .. math::

           E_{n_x,n_y,n_z}=\frac{\hbar^2\pi^2}{2m}
           \left(\frac{n_x^2}{L_x^2}+\frac{n_y^2}{L_y^2}
           +\frac{n_z^2}{L_z^2}\right).
        """
        if type(model) is not Piab3D:
            raise TypeError("model must be Piab3D")
        for name, value in (
            ("maximum_quantum_number_x", maximum_quantum_number_x),
            ("maximum_quantum_number_y", maximum_quantum_number_y),
            ("maximum_quantum_number_z", maximum_quantum_number_z),
        ):
            if type(value) is not int:
                raise TypeError(f"{name} must be a built-in integer")
            if value <= 0:
                raise ValueError(f"{name} must be positive")

        x_values = np.arange(1, maximum_quantum_number_x + 1, dtype=np.int64)
        y_values = np.arange(1, maximum_quantum_number_y + 1, dtype=np.int64)
        z_values = np.arange(1, maximum_quantum_number_z + 1, dtype=np.int64)
        n_x, n_y, n_z = np.meshgrid(
            x_values,
            y_values,
            z_values,
            indexing="ij",
        )
        triples = np.column_stack((n_x.ravel(), n_y.ravel(), n_z.ravel())).astype(
            np.int64
        )
        hbar = model.hbar_quantity
        mass = model.mass_quantity
        length_x = model.length_x_quantity
        length_y = model.length_y_quantity
        length_z = model.length_z_quantity
        energies = (
            hbar.magnitude**2
            * np.pi**2
            / (2.0 * mass.magnitude)
            * (
                triples[:, 0].astype(np.float64) ** 2 / length_x.magnitude**2
                + triples[:, 1].astype(np.float64) ** 2 / length_y.magnitude**2
                + triples[:, 2].astype(np.float64) ** 2 / length_z.magnitude**2
            )
        )
        if isinstance(model.unit_system.energy_unit, Unitless):
            energy_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(hbar.unit, PhysicalUnit):
                raise TypeError("physical model requires physical action units")
            if not isinstance(mass.unit, PhysicalUnit):
                raise TypeError("physical model requires physical mass units")
            if not isinstance(length_x.unit, PhysicalUnit):
                raise TypeError("physical model requires physical length units")
            source_unit = PhysicalUnit(
                expression=f"({hbar.unit.expression}) ** 2 / "
                f"({mass.unit.expression}) / "
                f"({length_x.unit.expression}) ** 2"
            )
            energy_unit = model.unit_system.energy_unit
            energies *= ModelSystemUnitConverter.conversion_factor(
                source_unit,
                energy_unit,
            )
        if not np.all(np.isfinite(energies)):
            raise ValueError("evaluated energies must be finite")
        order = np.lexsort(
            (
                triples[:, 2],
                triples[:, 1],
                triples[:, 0],
                energies,
            )
        )
        return Piab3DAnalyticalSolution(
            model=model,
            quantum_numbers=triples[order],
            energies=VectorQuantity(
                magnitude=energies[order],
                unit=energy_unit,
            ),
        )
