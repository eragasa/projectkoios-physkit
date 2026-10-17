"""Closed-form spectrum of the two-dimensional rectangular particle in a box."""

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

from ...base import Piab2D

type QuantumNumberPairs = npt.NDArray[np.int64]


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class Piab2DAnalyticalSolution(ResultsObject):
    """Retain one model and its ordered closed-form stationary-state energies.

    Parameters
    ----------
    model
        Exact rectangular particle-in-a-box model used for evaluation.
    quantum_numbers
        Unique positive integer pairs ``(n_x, n_y)`` with shape
        ``(state_count, 2)``. Pairs are ordered by nondecreasing energy, with
        lexicographic pair order breaking exact energy ties.
    energies
        One-dimensional energy quantity matching ``quantum_numbers`` and using
        the model's selected energy unit.

    Raises
    ------
    TypeError
        If a field has the wrong semantic type.
    ValueError
        If quantum-number shape, positivity, uniqueness, ordering, energy
        shape, or energy unit violates the result contract.
    """

    model: Piab2D
    quantum_numbers: QuantumNumberPairs
    energies: VectorQuantity

    def __post_init__(self) -> None:
        """Check each analytical-spectrum response argument."""
        self._check_arg_model()
        self._check_arg_quantum_numbers()
        self._check_arg_energies()

    def _check_arg_model(self) -> None:
        """Require the rectangular PIAB2D model that owns the spectrum."""
        if not isinstance(self.model, Piab2D):
            raise TypeError("model must be Piab2D")

    def _check_arg_quantum_numbers(self) -> None:
        """Require unique positive pairs in deterministic within-tie order."""
        if not isinstance(self.quantum_numbers, np.ndarray):
            raise TypeError("quantum_numbers must be a numpy.ndarray")
        values = np.asarray(self.quantum_numbers)
        if values.ndim != 2 or values.shape[1:] != (2,):
            raise ValueError("quantum_numbers must have shape (state_count, 2)")
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
        """Require model-unit eigenvalues ordered with their eigenstate pairs."""
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
class Piab2DAnalyticalEvaluator:
    """Evaluate a finite rectangular inventory of continuum energy levels."""

    def execute(
        self,
        *,
        model: Piab2D,
        maximum_quantum_number_x: int,
        maximum_quantum_number_y: int,
    ) -> Piab2DAnalyticalSolution:
        r"""Return ordered energies for selected positive quantum-number pairs.

        The returned inventory contains every pair satisfying
        ``1 <= n_x <= maximum_quantum_number_x`` and
        ``1 <= n_y <= maximum_quantum_number_y``. Energies are

        .. math::

           E_{n_x,n_y}=\frac{\hbar^2\pi^2}{2m}
           \left(\frac{n_x^2}{L_x^2}+\frac{n_y^2}{L_y^2}\right).

        Parameters
        ----------
        model
            Rectangular particle-in-a-box model.
        maximum_quantum_number_x
            Positive built-in integer upper bound for ``n_x``.
        maximum_quantum_number_y
            Positive built-in integer upper bound for ``n_y``.

        Returns
        -------
        Piab2DAnalyticalSolution
            Correlated model, ordered quantum-number pairs, and energies.

        Raises
        ------
        TypeError
            If an input has the wrong semantic type. Booleans and NumPy integer
            scalars are rejected for both bounds.
        ValueError
            If a bound is nonpositive or evaluated energies are nonfinite.
        """
        if type(model) is not Piab2D:
            raise TypeError("model must be Piab2D")
        for name, value in (
            ("maximum_quantum_number_x", maximum_quantum_number_x),
            ("maximum_quantum_number_y", maximum_quantum_number_y),
        ):
            if type(value) is not int:
                raise TypeError(f"{name} must be a built-in integer")
            if value <= 0:
                raise ValueError(f"{name} must be positive")

        x_values = np.arange(1, maximum_quantum_number_x + 1, dtype=np.int64)
        y_values = np.arange(1, maximum_quantum_number_y + 1, dtype=np.int64)
        n_x, n_y = np.meshgrid(x_values, y_values, indexing="ij")
        pairs = np.column_stack((n_x.ravel(), n_y.ravel())).astype(np.int64)
        hbar = model.hbar_quantity
        mass = model.mass_quantity
        length_x = model.length_x_quantity
        length_y = model.length_y_quantity
        energies = (
            hbar.magnitude**2
            * np.pi**2
            / (2.0 * mass.magnitude)
            * (
                pairs[:, 0].astype(np.float64) ** 2 / length_x.magnitude**2
                + pairs[:, 1].astype(np.float64) ** 2 / length_y.magnitude**2
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
        order = np.lexsort((pairs[:, 1], pairs[:, 0], energies))
        ordered_pairs = pairs[order]
        ordered_energies = energies[order]
        return Piab2DAnalyticalSolution(
            model=model,
            quantum_numbers=ordered_pairs,
            energies=VectorQuantity(
                magnitude=ordered_energies,
                unit=energy_unit,
            ),
        )
