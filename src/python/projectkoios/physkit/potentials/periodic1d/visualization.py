"""Visualization of one-dimensional periodic scalar potentials."""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from projectkoios.physkit.potentials.periodic.base import BasePeriodicPotential
from projectkoios.physkit.units import MatrixQuantity, Unitless, VectorQuantity


@dataclass(frozen=True, slots=True, kw_only=True)
class PeriodicPotential1DVisualizer:
    """Plot a one-dimensional periodic potential on fractional coordinates."""

    model: BasePeriodicPotential

    def __post_init__(self) -> None:
        """Require a one-dimensional periodic potential model."""
        if not isinstance(self.model, BasePeriodicPotential):
            raise TypeError("model must be BasePeriodicPotential")
        if self.model.dimension != 1:
            raise ValueError("model must be one-dimensional")

    def plot(
        self,
        *,
        fractional_coordinates: VectorQuantity,
        ax: Axes | None = None,
    ) -> Axes:
        """Plot potential values without displaying or closing the figure."""
        if not isinstance(fractional_coordinates, VectorQuantity):
            raise TypeError("fractional_coordinates must be VectorQuantity")
        if not isinstance(fractional_coordinates.unit, Unitless):
            raise ValueError("fractional_coordinates must be unitless")
        if fractional_coordinates.magnitude.size == 0:
            raise ValueError("fractional_coordinates must be nonempty")
        evaluation = self.model.evaluate(
            fractional_coordinates=MatrixQuantity(
                magnitude=fractional_coordinates.magnitude[:, np.newaxis],
                unit=Unitless(),
            )
        )
        if ax is None:
            _, ax = plt.subplots(figsize=(7.2, 3.8))
        ax.plot(
            fractional_coordinates.magnitude,
            evaluation.potential_values.magnitude,
        )
        ax.set_xlabel("Fractional coordinate $u$")
        ax.set_ylabel(f"Potential ({evaluation.potential_values.unit.expression})")
        ax.set_title("Periodic potential on one unit cell")
        ax.grid(alpha=0.25)
        return ax
