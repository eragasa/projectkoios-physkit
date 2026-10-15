# src/python/projectkoios/physkit/solidstate/electronic/results.py

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

type FloatArray = NDArray[np.float64]
type ComplexArray = NDArray[np.complex128]


@dataclass(frozen=True)
class ElectronicBandStructureResult:
    """
    Electronic band-structure solution.

    Parameters
    ----------
    k_points:
        Sampled Bloch wavevectors with shape ``(N_k,)``.
    energies:
        Band energies with shape ``(N_k, N_b)``.
    eigenvectors:
        Plane-wave coefficients with shape ``(N_k, N_G, N_b)``.
    reciprocal_vectors:
        Reciprocal basis vectors with shape ``(N_G,)``.
    """

    k_points: FloatArray
    energies: FloatArray
    eigenvectors: ComplexArray
    reciprocal_vectors: FloatArray

    @property
    def number_of_k_points(self) -> int:
        """
        Return the number of sampled wavevectors.

        Returns
        -------
        int
            Number of sampled wavevectors.
        """
        return self.k_points.size

    @property
    def number_of_bands(self) -> int:
        """
        Return the number of calculated bands.

        Returns
        -------
        int
            Number of calculated bands.
        """
        return self.energies.shape[1]

    def plane_wave_weights(
        self,
        k_index: int,
        band_index: int,
    ) -> FloatArray:
        """
        Return the plane-wave weights for one eigenstate.

        Parameters
        ----------
        k_index:
            Wavevector index.
        band_index:
            Band index.

        Returns
        -------
        FloatArray
            Weights ``|c_nk(G)|^2``.
        """
        coefficients = self.eigenvectors[
            k_index,
            :,
            band_index,
        ]

        return np.abs(coefficients) ** 2
