## Solver Types
# For TISE:
# - dense diagonalization;
# - sparse Lanczos;
# - Arnoldi;
# - shift-invert;
# - imaginary-time relaxation.

"""Deprecated shared TISE interfaces for one-dimensional PIAB solvers."""

from __future__ import annotations

from abc import ABC
from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from projectkoios.physkit.legacy.qm.models.piab1d import (
    Piab1D,
    Piab1DResult,
    Piab1DSolver,
)

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


class Piab1DTISEResult(
    Piab1DResult,
    ABC,
):
    """
    Base result produced by a Piab1D TISE solver.

    Parameters
    ----------
    model:
        Particle-in-a-box model solved by the solver.
    quantum_numbers:
        Positive integer quantum numbers associated with the stationary
        states.
    energies:
        Energy eigenvalues ordered by increasing quantum number.

    Attributes
    ----------
    model:
        Physical model associated with the result.
    quantum_numbers:
        Positive integer quantum numbers.
    energies:
        Ordered energy eigenvalues.
    number_of_states:
        Number of stationary states contained in the result.

    Notes
    -----
    This class defines the information shared by analytical and numerical
    TISE results. The representation of the stationary eigenfunctions is
    defined by concrete subclasses.
    """

    def __init__(
        self,
        model: Piab1D,
        quantum_numbers: IntArray,
        energies: FloatArray,
    ) -> None:
        if not isinstance(model, Piab1D):
            raise TypeError("model must be an instance of Piab1D.")

        quantum_numbers_array = np.asarray(
            quantum_numbers,
        )

        energies_array = np.asarray(
            energies,
            dtype=np.float64,
        )

        self.check_arrays(
            quantum_numbers=quantum_numbers_array,
            energies=energies_array,
        )

        self.model: Piab1D = model

        self.quantum_numbers: IntArray = np.array(
            quantum_numbers_array,
            dtype=np.int64,
            copy=True,
        )

        self.energies: FloatArray = np.array(
            energies_array,
            dtype=np.float64,
            copy=True,
        )

        self.quantum_numbers.setflags(
            write=False,
        )

        self.energies.setflags(
            write=False,
        )

    @staticmethod
    def check_arrays(
        quantum_numbers: NDArray[np.generic],
        energies: FloatArray,
    ) -> None:
        """
        Validate the quantum numbers and energy eigenvalues.

        Parameters
        ----------
        quantum_numbers:
            Candidate quantum-number array.
        energies:
            Candidate energy array.

        Raises
        ------
        TypeError
            If the quantum numbers do not have an integer dtype.
        ValueError
            If either array has an invalid shape, the array lengths differ,
            the quantum numbers are not consecutive positive integers, or
            the energies are nonfinite or not strictly increasing.
        """
        if quantum_numbers.ndim != 1:
            raise ValueError("quantum_numbers must be one-dimensional.")

        if energies.ndim != 1:
            raise ValueError("energies must be one-dimensional.")

        if quantum_numbers.size == 0:
            raise ValueError("At least one quantum number is required.")

        if quantum_numbers.size != energies.size:
            raise ValueError("quantum_numbers and energies must have the same length.")

        if not np.issubdtype(
            quantum_numbers.dtype,
            np.integer,
        ):
            raise TypeError("quantum_numbers must have an integer dtype.")

        expected_quantum_numbers = np.arange(
            1,
            quantum_numbers.size + 1,
            dtype=np.int64,
        )

        if not np.array_equal(
            quantum_numbers,
            expected_quantum_numbers,
        ):
            raise ValueError(
                "quantum_numbers must be consecutive positive "
                "integers beginning with one."
            )

        if not np.all(np.isfinite(energies)):
            raise ValueError("energies must contain only finite values.")

        if np.any(energies <= 0.0):
            raise ValueError("energies must be greater than zero.")

        if energies.size > 1 and np.any(np.diff(energies) <= 0.0):
            raise ValueError("energies must be strictly increasing.")

    @property
    def number_of_states(self) -> int:
        """
        Return the number of stationary states.

        Returns
        -------
        int
            Number of states contained in the result.
        """
        return int(self.energies.size)

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}("
            f"model={self.model!r}, "
            f"number_of_states={self.number_of_states})"
        )


class Piab1DTISESolver[Piab1DTISEResultType: Piab1DTISEResult](
    Piab1DSolver[Piab1DTISEResultType],
    ABC,
):
    """
    Base interface for a Piab1D TISE solver.

    Parameters
    ----------
    model:
        Physical particle-in-a-box model.
    number_of_states:
        Number of lowest-energy stationary states to compute.

    Attributes
    ----------
    model:
        Physical model supplied to the solver.
    number_of_states:
        Number of stationary states requested.

    Notes
    -----
    Construction stores and validates the solver configuration. Computation
    is performed only when ``solve()`` is called.
    """

    def __init__(
        self,
        model: Piab1D,
        number_of_states: int,
    ) -> None:
        if not isinstance(model, Piab1D):
            raise TypeError("model must be an instance of Piab1D.")

        if isinstance(
            number_of_states,
            (bool, np.bool_),
        ) or not isinstance(
            number_of_states,
            Integral,
        ):
            raise TypeError("number_of_states must be an integer.")

        if number_of_states < 1:
            raise ValueError("number_of_states must be at least one.")

        self.model: Piab1D = model
        self.number_of_states: int = int(number_of_states)
