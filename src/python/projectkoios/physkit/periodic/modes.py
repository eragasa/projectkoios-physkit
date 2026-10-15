# src/python/projectkoios/physkit/periodic/modes.py

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from projectkoios.physkit.periodic.basis import ReciprocalModeBasis1D

type FloatArray = NDArray[np.float64]
type ComplexArray = NDArray[np.complex128]


class BlochMode1D:
    """
    Bloch mode represented in a reciprocal-lattice Fourier basis.

    A Bloch mode has the form

    .. math::

        f_k(x)
        =
        e^{ikx}u_k(x),

    where the cell-periodic part is

    .. math::

        u_k(x)
        =
        \\sum_G c_k(G)e^{iGx}.

    The object represents the mathematical structure of a Bloch mode. It
    does not assign a physical interpretation to the coefficients or field.
    The mode may describe an electronic wavefunction, a lattice-displacement
    field, an electromagnetic field, or another periodic quantity.

    Parameters
    ----------
    wavevector:
        Bloch wavevector ``k``.
    basis:
        Reciprocal-mode basis defining the vectors ``G``.
    coefficients:
        Complex Fourier coefficients ``c_k(G)``. The array must contain one
        coefficient for every reciprocal mode in ``basis``.

    Attributes
    ----------
    wavevector:
        Bloch wavevector.
    basis:
        Reciprocal-mode basis.
    coefficients:
        Read-only complex Fourier coefficients with shape
        ``(number_of_modes,)``.
    """

    def __init__(
        self,
        wavevector: float,
        basis: ReciprocalModeBasis1D,
        coefficients: ComplexArray,
    ) -> None:
        self.wavevector: float = float(wavevector)

        self.basis = basis

        # Normalize the coefficient representation so that all subsequent
        # mode reconstruction uses a consistent complex dtype.
        self.coefficients: ComplexArray = np.asarray(
            coefficients,
            dtype=np.complex128,
        )

        self.check_args()

        # A Bloch mode represents one fixed state. Preventing in-place changes
        # keeps the basis and coefficient vector synchronized.
        self.coefficients.setflags(write=False)

    def check_args(self) -> None:
        """
        Validate the Bloch-mode arguments.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If the wavevector is not finite, the coefficient array is not
            one-dimensional, or its size does not match the basis.
        """
        if not np.isfinite(self.wavevector):
            raise ValueError("wavevector must be finite.")

        if self.coefficients.ndim != 1:
            raise ValueError("coefficients must be one-dimensional.")

        if self.coefficients.size != self.basis.size:
            raise ValueError(
                "The number of coefficients must equal the reciprocal-mode basis size."
            )

    def cell_periodic_part(
        self,
        coordinates: FloatArray,
    ) -> ComplexArray:
        """
        Reconstruct the cell-periodic part of the Bloch mode.

        The reconstruction is

        .. math::

            u_k(x)
            =
            \\sum_G c_k(G)e^{iGx}.

        Parameters
        ----------
        coordinates:
            One-dimensional real-space coordinates with shape
            ``(number_of_coordinates,)``.

        Returns
        -------
        ComplexArray
            Values of ``u_k(x)`` at the supplied coordinates.
        """
        coordinates = np.asarray(
            coordinates,
            dtype=np.float64,
        )

        if coordinates.ndim != 1:
            raise ValueError("coordinates must be one-dimensional.")

        # Each row contains one reciprocal-lattice basis function evaluated
        # over every real-space coordinate:
        #
        #     plane_waves[m, j] = exp(i G_m x_j).
        plane_waves: ComplexArray = np.exp(
            1j * self.basis.vectors[:, None] * coordinates[None, :]
        )

        # Contract the reciprocal-mode index to obtain u_k(x_j).
        return self.coefficients @ plane_waves

    def values(
        self,
        coordinates: FloatArray,
    ) -> ComplexArray:
        """
        Reconstruct the complete Bloch mode.

        The complete mode is

        .. math::

            f_k(x)
            =
            e^{ikx}u_k(x).

        Parameters
        ----------
        coordinates:
            One-dimensional real-space coordinates.

        Returns
        -------
        ComplexArray
            Values of the complete Bloch mode.
        """
        coordinates = np.asarray(
            coordinates,
            dtype=np.float64,
        )

        periodic_part = self.cell_periodic_part(coordinates)

        # Multiply the cell-periodic function by the macroscopic Bloch phase.
        bloch_phase: ComplexArray = np.exp(1j * self.wavevector * coordinates)

        return bloch_phase * periodic_part

    def reciprocal_weights(
        self,
    ) -> FloatArray:
        """
        Calculate the normalized reciprocal-mode weights.

        Returns
        -------
        FloatArray
            Values ``|c_k(G)|^2`` for every reciprocal mode.

        Notes
        -----
        These weights describe the representation of the generic Bloch mode
        in the reciprocal basis. Their physical interpretation depends on
        the domain that produced the mode.
        """
        return np.abs(self.coefficients) ** 2
