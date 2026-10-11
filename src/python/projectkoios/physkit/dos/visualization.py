"""Dimension-independent visualization of density-of-states results."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.dos.base import (
    GaussianBandDensityOfStates,
    HistogramBandDensityOfStates,
)
from projectkoios.physkit.dos.thermal import FermiDiracThermalDensityOfStates
from projectkoios.physkit.units import VectorQuantity

type DensityOfStatesOrientation = Literal["energy_horizontal", "energy_vertical"]
type ThermalDensityOfStatesComponent = Literal[
    "occupied",
    "unoccupied",
    "thermally_active",
]


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class DensityOfStatesPlotData(DataObject):
    """Represent correlated energy coordinates and DOS values."""

    energies: VectorQuantity
    density_of_states: VectorQuantity

    def __post_init__(self) -> None:
        """Check nonempty correlated plotting vectors."""
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if not isinstance(self.density_of_states, VectorQuantity):
            raise TypeError("density_of_states must be VectorQuantity")
        if self.energies.magnitude.size == 0:
            raise ValueError("energies must be nonempty")
        if self.density_of_states.magnitude.shape != self.energies.magnitude.shape:
            raise ValueError("density_of_states must match energies")
        if np.any(self.density_of_states.magnitude < 0.0):
            raise ValueError("density_of_states must be nonnegative")

    @classmethod
    def from_gaussian(
        cls,
        result: GaussianBandDensityOfStates,
    ) -> DensityOfStatesPlotData:
        """Adapt a Gaussian DOS result without copying its quantity records."""
        if not isinstance(result, GaussianBandDensityOfStates):
            raise TypeError("result must be GaussianBandDensityOfStates")
        return cls(
            energies=result.request.energy_grid,
            density_of_states=result.density_of_states,
        )

    @classmethod
    def from_histogram(
        cls,
        result: HistogramBandDensityOfStates,
    ) -> DensityOfStatesPlotData:
        """Adapt a histogram DOS result without copying its quantity records."""
        if not isinstance(result, HistogramBandDensityOfStates):
            raise TypeError("result must be HistogramBandDensityOfStates")
        return cls(
            energies=result.energy_bin_centers,
            density_of_states=result.density_of_states,
        )

    @classmethod
    def from_thermal(
        cls,
        result: FermiDiracThermalDensityOfStates,
        *,
        component: ThermalDensityOfStatesComponent,
    ) -> DensityOfStatesPlotData:
        """Adapt one occupation-resolved thermal DOS component.

        Parameters
        ----------
        result
            Correlated thermal-DOS result.
        component
            Occupied, unoccupied, or thermally active density to render.

        Returns
        -------
        DensityOfStatesPlotData
            Energy coordinates paired with the selected density vector.
        """
        if not isinstance(result, FermiDiracThermalDensityOfStates):
            raise TypeError("result must be FermiDiracThermalDensityOfStates")
        if component == "occupied":
            density = result.occupied_density_of_states
        elif component == "unoccupied":
            density = result.unoccupied_density_of_states
        elif component == "thermally_active":
            density = result.thermally_active_density_of_states
        else:
            raise ValueError("component must identify a supported thermal DOS")
        return cls(
            energies=result.request.energies,
            density_of_states=density,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class DensityOfStatesVisualizer:
    """Plot a DOS without assuming the reciprocal-space dimension."""

    data: DensityOfStatesPlotData

    def __post_init__(self) -> None:
        """Require complete DOS plot data."""
        if not isinstance(self.data, DensityOfStatesPlotData):
            raise TypeError("data must be DensityOfStatesPlotData")

    def plot(
        self,
        *,
        ax: Axes | None = None,
        orientation: DensityOfStatesOrientation = "energy_horizontal",
        label: str | None = None,
        title: str = "Density of states",
    ) -> Axes:
        """Plot energy horizontally or vertically and return the axes."""
        if ax is not None and not isinstance(ax, Axes):
            raise TypeError("ax must be matplotlib.axes.Axes or None")
        if orientation not in ("energy_horizontal", "energy_vertical"):
            raise ValueError("orientation must identify a supported DOS layout")
        if label is not None and not isinstance(label, str):
            raise TypeError("label must be a string or None")
        if not isinstance(title, str):
            raise TypeError("title must be a string")
        if ax is None:
            _, ax = plt.subplots(figsize=(6.8, 4.3))
        energies = self.data.energies.magnitude
        density = self.data.density_of_states.magnitude
        energy_label = f"Energy ({self.data.energies.unit.expression})"
        density_label = (
            "Density of states "
            f"({self.data.density_of_states.unit.expression})"
        )
        if orientation == "energy_horizontal":
            ax.plot(energies, density, label=label)
            ax.set_xlabel(energy_label)
            ax.set_ylabel(density_label)
        else:
            ax.plot(density, energies, label=label)
            ax.set_xlabel(density_label)
            ax.set_ylabel(energy_label)
        ax.set_title(title)
        ax.grid(alpha=0.25)
        return ax
