# src/python/projectkoios/physkit/solidstate/visualization/phonon_bands.py

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from projectkoios.physkit.periodic.reciprocal import ReciprocalLattice1D
from projectkoios.physkit.solidstate.phonons.results import (
    PhononBandStructureResult,
)


class PhononBandStructureVisualizer:
    """
    Visualizer for a one-dimensional phonon dispersion.
    """

    def __init__(
        self,
        result: PhononBandStructureResult,
        reciprocal_lattice: ReciprocalLattice1D,
    ) -> None:
        self.result = result
        self.reciprocal_lattice = reciprocal_lattice

    def plot(
        self,
        ax: Axes | None = None,
    ) -> Axes:
        """
        Plot the phonon angular-frequency branches.

        Parameters
        ----------
        ax:
            Existing Matplotlib axes. If omitted, new axes are created.

        Returns
        -------
        Axes
            Axes containing the phonon dispersion.
        """
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))

        normalized_q_points = (
            self.result.q_points / self.reciprocal_lattice.fundamental_vector
        )

        for branch_index in range(self.result.number_of_branches):
            ax.plot(
                normalized_q_points,
                self.result.angular_frequencies[
                    :,
                    branch_index,
                ],
                linewidth=2.0,
            )

        ax.set_xlabel(r"Phonon wavevector $q/G_0$")

        ax.set_ylabel(r"Angular frequency $\omega_\nu(q)$")

        ax.set_title("Phonon Dispersion")

        ax.grid(alpha=0.25)

        return ax
