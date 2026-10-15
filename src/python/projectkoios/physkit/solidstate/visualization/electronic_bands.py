# src/python/projectkoios/physkit/solidstate/visualization/electronic_bands.py

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from projectkoios.physkit.periodic.reciprocal import ReciprocalLattice1D
from projectkoios.physkit.solidstate.electronic.results import (
    ElectronicBandStructureResult,
)


class ElectronicBandStructureVisualizer:
    """
    Visualizer for a one-dimensional electronic band structure.
    """

    def __init__(
        self,
        result: ElectronicBandStructureResult,
        reciprocal_lattice: ReciprocalLattice1D,
    ) -> None:
        self.result = result
        self.reciprocal_lattice = reciprocal_lattice

    def plot(
        self,
        ax: Axes | None = None,
    ) -> Axes:
        """
        Plot the calculated electronic bands.

        Parameters
        ----------
        ax:
            Existing Matplotlib axes. If omitted, new axes are created.

        Returns
        -------
        Axes
            Axes containing the band-structure plot.
        """
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))

        normalized_k_points = (
            self.result.k_points / self.reciprocal_lattice.fundamental_vector
        )

        for band_index in range(self.result.number_of_bands):
            ax.plot(
                normalized_k_points,
                self.result.energies[
                    :,
                    band_index,
                ],
                linewidth=2.0,
            )

        ax.axvline(
            -0.5,
            color="black",
            linestyle="--",
            linewidth=1.0,
        )

        ax.axvline(
            0.5,
            color="black",
            linestyle="--",
            linewidth=1.0,
        )

        ax.set_xlabel(r"Bloch wavevector $k/G_0$")

        ax.set_ylabel(r"Energy $E_n(k)$")

        ax.set_title("Electronic Band Structure")

        ax.grid(alpha=0.25)

        return ax
