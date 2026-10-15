"""Unit-aware interval representations migrated from ksdft2effmass."""

from __future__ import annotations

from dataclasses import dataclass

from .boundary_conditions import DirichletBoundaryCondition
from .cartesian_grids import UniformCartesianGrid1D


@dataclass(frozen=True, slots=True)
class DirichletInterval:
    """Compose a uniform one-dimensional grid and Dirichlet boundary data."""

    grid: UniformCartesianGrid1D
    boundary_condition: DirichletBoundaryCondition

    def __post_init__(self) -> None:
        if not isinstance(self.grid, UniformCartesianGrid1D):
            raise TypeError("grid must be UniformCartesianGrid1D")
        if not isinstance(
            self.boundary_condition,
            DirichletBoundaryCondition,
        ):
            raise TypeError("boundary_condition must be DirichletBoundaryCondition")
        if self.grid.interior_point_count <= 0:
            raise ValueError("Dirichlet interval must contain an interior point")

    @property
    def interior_points(self) -> int:
        """Return the number of ordered interior coordinates."""
        return self.grid.interior_point_count
