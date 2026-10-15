# src/python/projectkoios/physkit/solidstate/phonons/solvers.py

from __future__ import annotations

import numpy as np

from projectkoios.physkit.periodic.kpoints import KPointPath1D
from projectkoios.physkit.solidstate.phonons.operators import (
    MonatomicDynamicalMatrix1D,
)
from projectkoios.physkit.solidstate.phonons.results import (
    PhononBandStructureResult,
)


class PhononBlochSolver1D:
    """
    Bloch eigensolver for a one-dimensional phonon model.

    Parameters
    ----------
    dynamical_matrix:
        Dynamical-matrix builder for the phonon model.
    """

    def __init__(
        self,
        dynamical_matrix: MonatomicDynamicalMatrix1D,
    ) -> None:
        self.dynamical_matrix = dynamical_matrix

    def solve(
        self,
        q_path: KPointPath1D,
    ) -> PhononBandStructureResult:
        """
        Calculate the phonon dispersion over a wavevector path.

        Parameters
        ----------
        q_path:
            Ordered path of phonon wavevectors.

        Returns
        -------
        PhononBandStructureResult
            Frequencies and polarization vectors.
        """
        number_of_q_points = q_path.size

        angular_frequencies = np.zeros(
            (number_of_q_points, 1),
            dtype=np.float64,
        )

        polarization_vectors = np.zeros(
            (number_of_q_points, 1, 1),
            dtype=np.complex128,
        )

        for q_index, phonon_wavevector in enumerate(q_path.values):
            matrix = self.dynamical_matrix.construct_matrix(phonon_wavevector)

            squared_frequencies, eigenvectors = np.linalg.eigh(matrix)

            # Small negative eigenvalues may occur through floating-point
            # roundoff at an exact acoustic zero mode.
            squared_frequencies = np.maximum(
                squared_frequencies,
                0.0,
            )

            angular_frequencies[q_index] = np.sqrt(squared_frequencies)

            polarization_vectors[q_index] = eigenvectors

        return PhononBandStructureResult(
            q_points=q_path.values.copy(),
            angular_frequencies=angular_frequencies,
            polarization_vectors=polarization_vectors,
        )
