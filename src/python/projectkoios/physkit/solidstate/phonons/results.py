# src/python/projectkoios/physkit/solidstate/phonons/results.py

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

type FloatArray = NDArray[np.float64]
type ComplexArray = NDArray[np.complex128]


@dataclass(frozen=True)
class PhononBandStructureResult:
    """
    Phonon dispersion solution.

    Parameters
    ----------
    q_points:
        Sampled phonon wavevectors with shape ``(N_q,)``.
    angular_frequencies:
        Phonon angular frequencies with shape ``(N_q, N_b)``.
    polarization_vectors:
        Phonon eigenvectors with shape ``(N_q, N_d, N_b)``.
    """

    q_points: FloatArray
    angular_frequencies: FloatArray
    polarization_vectors: ComplexArray

    @property
    def number_of_q_points(self) -> int:
        """
        Return the number of sampled phonon wavevectors.

        Returns
        -------
        int
            Number of sampled wavevectors.
        """
        return self.q_points.size

    @property
    def number_of_branches(self) -> int:
        """
        Return the number of phonon branches.

        Returns
        -------
        int
            Number of phonon branches.
        """
        return self.angular_frequencies.shape[1]
