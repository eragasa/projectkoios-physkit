"""Dimension-independent visualization of sampled band structures."""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.units import (
    MatrixQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    Unitless,
    VectorQuantity,
)

_ENERGY_UNIT = PhysicalUnit("joule")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class BandStructurePlotData(DataObject):
    """Represent band energies along an ordered scalar path coordinate."""

    path_coordinates: VectorQuantity
    band_energies: MatrixQuantity

    def __post_init__(self) -> None:
        """Check path and band-energy correlation."""
        if not isinstance(self.path_coordinates, VectorQuantity):
            raise TypeError("path_coordinates must be VectorQuantity")
        if not isinstance(self.band_energies, MatrixQuantity):
            raise TypeError("band_energies must be MatrixQuantity")
        if self.path_coordinates.magnitude.size == 0:
            raise ValueError("path_coordinates must be nonempty")
        if isinstance(self.band_energies.unit, PhysicalUnit) and not (
            ModelSystemUnitConverter.compatible(
                self.band_energies.unit,
                _ENERGY_UNIT,
            )
        ):
            raise ValueError("band_energies must use units compatible with energy")
        if not isinstance(self.band_energies.unit, PhysicalUnit | Unitless):
            raise TypeError("band_energies must use PhysicalUnit or Unitless")
        if self.band_energies.magnitude.shape[0] != (
            self.path_coordinates.magnitude.size
        ):
            raise ValueError("band_energies rows must match path coordinates")
        if self.band_energies.magnitude.shape[1] == 0:
            raise ValueError("band_energies must contain at least one band")


@dataclass(frozen=True, slots=True, kw_only=True)
class BandStructureVisualizer:
    """Plot bands along a caller-defined path without assuming dimension."""

    data: BandStructurePlotData

    def __post_init__(self) -> None:
        """Require a complete dimension-independent plot-data record."""
        if not isinstance(self.data, BandStructurePlotData):
            raise TypeError("data must be BandStructurePlotData")

    def plot(
        self,
        *,
        ax: Axes | None = None,
        reference_band_energies: MatrixQuantity | None = None,
        path_label: str = "Band-path coordinate",
        title: str = "Band structure",
    ) -> Axes:
        """Plot primary bands and an optional same-grid reference."""
        if ax is not None and not isinstance(ax, Axes):
            raise TypeError("ax must be matplotlib.axes.Axes or None")
        if not isinstance(path_label, str):
            raise TypeError("path_label must be a string")
        if not isinstance(title, str):
            raise TypeError("title must be a string")
        if reference_band_energies is not None:
            self._check_reference_band_energies(reference_band_energies)
        if ax is None:
            _, ax = plt.subplots(figsize=(7.4, 4.8))
        coordinates = self.data.path_coordinates.magnitude
        if reference_band_energies is not None:
            ax.plot(
                coordinates,
                reference_band_energies.magnitude,
                color="0.65",
                linestyle="--",
                linewidth=0.9,
            )
        ax.plot(
            coordinates,
            self.data.band_energies.magnitude,
            linewidth=1.2,
        )
        ax.set_xlabel(path_label)
        ax.set_ylabel(f"Energy ({self.data.band_energies.unit.expression})")
        ax.set_title(title)
        ax.grid(alpha=0.25)
        return ax

    def _check_reference_band_energies(
        self,
        reference_band_energies: MatrixQuantity,
    ) -> None:
        """Require reference bands in the primary energy representation."""
        if not isinstance(reference_band_energies, MatrixQuantity):
            raise TypeError("reference_band_energies must be MatrixQuantity or None")
        if reference_band_energies.unit != self.data.band_energies.unit:
            raise ValueError("reference_band_energies must use the primary energy unit")
        if (
            reference_band_energies.magnitude.shape
            != self.data.band_energies.magnitude.shape
        ):
            raise ValueError(
                "reference_band_energies must match the primary band shape"
            )
