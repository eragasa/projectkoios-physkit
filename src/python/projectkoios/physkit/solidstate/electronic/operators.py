# src/python/projectkoios/physkit/solidstate/electronic/operators.py

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from projectkoios.physkit.periodic.basis import ReciprocalModeBasis1D
from projectkoios.physkit.solidstate.electronic.models import (
    PlaneWaveBlochModel1D,
)

type FloatArray = NDArray[np.float64]
type ComplexArray = NDArray[np.complex128]


class PlaneWaveBlochHamiltonian1D:
    """
    Electronic Bloch Hamiltonian in a reciprocal-mode basis.

    The matrix elements are

    .. math::

        H_{G'G}(k)
        =
        \\frac{\\hbar^2}{2m}(k+G)^2\\delta_{G'G}
        +
        V_{G'-G}.
    """

    def __init__(
        self,
        model: PlaneWaveBlochModel1D,
        basis: ReciprocalModeBasis1D,
    ) -> None:
        self.model = model
        self.basis = basis

        # The potential does not depend on k, so it can be constructed once
        # and reused for every wavevector in the Brillouin zone.
        self.potential_matrix: ComplexArray = self.construct_potential_matrix()

    def kinetic_spectrum(
        self,
        bloch_wavevector: float,
    ) -> FloatArray:
        """
        Calculate the kinetic-energy eigenvalues.

        Parameters
        ----------
        bloch_wavevector:
            Electronic Bloch wavevector ``k``.

        Returns
        -------
        FloatArray
            Diagonal kinetic-energy entries.
        """
        shifted_wavevectors = self.basis.shifted_wavevectors(bloch_wavevector)

        return self.model.hbar**2 * shifted_wavevectors**2 / (2.0 * self.model.mass)

    def construct_potential_matrix(
        self,
    ) -> ComplexArray:
        """
        Construct the plane-wave potential matrix.

        Returns
        -------
        ComplexArray
            Matrix with elements ``V_(G' - G)``.
        """
        transferred_wavevectors = (
            self.basis.vectors[:, None] - self.basis.vectors[None, :]
        )

        return self.model.potential.coefficients(transferred_wavevectors)

    def construct_matrix(
        self,
        bloch_wavevector: float,
    ) -> ComplexArray:
        """
        Construct the Bloch Hamiltonian at one wavevector.

        Parameters
        ----------
        bloch_wavevector:
            Bloch wavevector ``k``.

        Returns
        -------
        ComplexArray
            Hermitian Hamiltonian matrix ``H(k)``.
        """
        hamiltonian = self.potential_matrix.copy()

        # Add the kinetic spectrum directly to the diagonal without
        # constructing a separate dense diagonal matrix.
        diagonal_indices = np.diag_indices(self.basis.size)

        hamiltonian[diagonal_indices] += self.kinetic_spectrum(bloch_wavevector)

        return hamiltonian
