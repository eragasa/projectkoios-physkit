"""Deprecated analytical TISE solution of the one-dimensional PIAB."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from projectkoios.physkit.legacy.qm.solvers.piab1d.tise import (
    Piab1DTISEResult,
    Piab1DTISESolver,
)

FloatArray = NDArray[np.float64]


class Piab1DAnalyticalTISEResult(
    Piab1DTISEResult,
):
    """
    Analytical TISE result for a one-dimensional particle in a box.

    The energy eigenvalues are

    .. math::

        E_n
        =
        \\frac{
            \\hbar^2 \\pi^2 n^2
        }{
            2mL^2
        },

    where

    .. math::

        L
        =
        x_{\\mathrm{upper}}
        -
        x_{\\mathrm{lower}}.

    The normalized stationary eigenfunctions are

    .. math::

        \\psi_n(x)
        =
        \\sqrt{\\frac{2}{L}}
        \\sin\\left[
            \\frac{
                n\\pi
                (x-x_{\\mathrm{lower}})
            }{
                L
            }
        \\right].

    Notes
    -----
    The result stores the model, quantum numbers, and energies through
    ``Piab1DTISEResult``. Eigenfunctions are evaluated when requested
    because the analytical result is not tied to a numerical grid.
    """

    def evaluate_eigenfunctions(
        self,
        coordinates: FloatArray,
    ) -> FloatArray:
        """
        Evaluate all stored stationary eigenfunctions.

        Parameters
        ----------
        coordinates:
            One-dimensional coordinate array. Every coordinate must lie
            within the closed interval of the model.

        Returns
        -------
        numpy.ndarray
            Eigenfunction values with shape

            ``(number_of_coordinates, number_of_states)``.

            Column ``state_index`` contains the eigenfunction associated
            with ``quantum_numbers[state_index]``.

        Raises
        ------
        TypeError
            If ``coordinates`` is complex-valued.
        ValueError
            If ``coordinates`` is not one-dimensional, contains nonfinite
            values, or extends outside the model domain.
        """
        coordinates_input = np.asarray(
            coordinates,
        )

        if not np.isrealobj(coordinates_input):
            raise TypeError("coordinates must be real-valued.")

        coordinates_array = np.asarray(
            coordinates_input,
            dtype=np.float64,
        )

        if coordinates_array.ndim != 1:
            raise ValueError("coordinates must be one-dimensional.")

        if not np.all(np.isfinite(coordinates_array)):
            raise ValueError("coordinates must contain only finite values.")

        if np.any(coordinates_array < self.model.x_lower):
            raise ValueError("coordinates must not be less than model.x_lower.")

        if np.any(coordinates_array > self.model.x_upper):
            raise ValueError("coordinates must not be greater than model.x_upper.")

        relative_coordinates = coordinates_array - self.model.x_lower

        wave_numbers = np.pi * self.quantum_numbers / self.model.length

        eigenfunctions = np.sqrt(2.0 / self.model.length) * np.sin(
            relative_coordinates[:, np.newaxis] * wave_numbers[np.newaxis, :]
        )

        result = np.asarray(
            eigenfunctions,
            dtype=np.float64,
        )

        result.setflags(
            write=False,
        )

        return result

    def evaluate_eigenfunction(
        self,
        coordinates: FloatArray,
        quantum_number: int,
    ) -> FloatArray:
        """
        Evaluate one stationary eigenfunction.

        Parameters
        ----------
        coordinates:
            One-dimensional coordinate array.
        quantum_number:
            Positive integer quantum number contained in this result.

        Returns
        -------
        numpy.ndarray
            Eigenfunction values with shape
            ``(number_of_coordinates,)``.

        Raises
        ------
        TypeError
            If ``quantum_number`` is not an integer.
        ValueError
            If the requested quantum number is not contained in the
            result.
        """
        if isinstance(
            quantum_number,
            (bool, np.bool_),
        ) or not isinstance(
            quantum_number,
            (int, np.integer),
        ):
            raise TypeError("quantum_number must be an integer.")

        matching_indices = np.flatnonzero(self.quantum_numbers == quantum_number)

        if matching_indices.size == 0:
            raise ValueError(
                f"quantum_number={quantum_number} is not contained in this result."
            )

        state_index = int(matching_indices[0])

        eigenfunctions = self.evaluate_eigenfunctions(
            coordinates=coordinates,
        )

        eigenfunction = np.array(
            eigenfunctions[:, state_index],
            dtype=np.float64,
            copy=True,
        )

        eigenfunction.setflags(
            write=False,
        )

        return eigenfunction


class Piab1DAnalyticalTISESolver(
    Piab1DTISESolver[Piab1DAnalyticalTISEResult],
):
    """
    Analytical TISE solver for a one-dimensional particle in a box.

    Construction stores the model and number of requested states.
    Computation occurs only when ``solve()`` is called.
    """

    def solve(
        self,
    ) -> Piab1DAnalyticalTISEResult:
        """
        Compute the analytical stationary-state energies.

        Returns
        -------
        Piab1DAnalyticalTISEResult
            Analytical energy spectrum and eigenfunction evaluator.
        """
        quantum_numbers = np.arange(
            1,
            self.number_of_states + 1,
            dtype=np.int64,
        )

        hbar = self.model.constants.hbar
        mass = self.model.mass
        length = self.model.length

        energies = (
            hbar**2
            * np.pi**2
            * quantum_numbers.astype(np.float64) ** 2
            / (2.0 * mass * length**2)
        )

        return Piab1DAnalyticalTISEResult(
            model=self.model,
            quantum_numbers=quantum_numbers,
            energies=np.asarray(
                energies,
                dtype=np.float64,
            ),
        )
