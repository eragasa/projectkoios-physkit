# src/python/projectkoios/physkit/solidstate/visualization/bloch_modes.py

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from projectkoios.physkit.solidstate.electronic.results import (
    ElectronicBandStructureResult,
)


class ElectronicBlochModeVisualizer1D:
    """
    Visualizer for electronic Bloch eigenstates.
    """

    def __init__(
        self,
        result: ElectronicBandStructureResult,
    ) -> None:
        self.result = result

    def plot_plane_wave_weights(
        self,
        k_index: int,
        band_index: int,
        ax: Axes | None = None,
    ) -> Axes:
        """
        Plot the reciprocal-space composition of a Bloch state.

        Parameters
        ----------
        k_index:
            Index of the selected Bloch wavevector.
        band_index:
            Index of the selected electronic band.
        ax:
            Existing Matplotlib axes.

        Returns
        -------
        Axes
            Axes containing the plane-wave weights.
        """
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 5))

        weights = self.result.plane_wave_weights(
            k_index,
            band_index,
        )

        ax.bar(
            self.result.reciprocal_vectors,
            weights,
            width=0.8,
            color="tab:blue",
            edgecolor="black",
        )

        ax.set_xlabel(r"Reciprocal vector $G$")

        ax.set_ylabel(r"Weight $|c_{nk}(G)|^2$")

        ax.set_title("Plane-Wave Composition")

        ax.grid(
            axis="y",
            alpha=0.25,
        )

        return ax

    def reconstruct_cell_periodic_state(
        self,
        x,
        *,
        k_index: int,
        band_index: int,
    ):
        """
        Reconstruct the cell-periodic Bloch function.

        Parameters
        ----------
        x:
            Real-space coordinates.
        k_index:
            Index of the selected Bloch wavevector.
        band_index:
            Index of the selected band.

        Returns
        -------
        NDArray
            Complex values of ``u_nk(x)``.
        """
        coordinates = np.asarray(
            x,
            dtype=np.float64,
        )

        coefficients = self.result.eigenvectors[
            k_index,
            :,
            band_index,
        ]

        plane_waves = np.exp(
            1j * self.result.reciprocal_vectors[:, None] * coordinates[None, :]
        )

        return coefficients @ plane_waves
