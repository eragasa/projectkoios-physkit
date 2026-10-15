"""Unit-aware uniform Cartesian grids migrated from ksdft2effmass."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    ModelSystemUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True)
class UniformCartesianGrid1D:
    """Represent one ordered uniform Cartesian axis including both boundaries.

    Parameters
    ----------
    lower_bound:
        Finite lower coordinate with an explicit physical or unitless unit.
    upper_bound:
        Finite upper coordinate compatible with ``lower_bound``.
    requested_spacing:
        Positive compatible spacing. The domain length must be an integer
        multiple of this spacing to absolute tolerance ``1e-12``.
    """

    lower_bound: ScalarQuantity
    upper_bound: ScalarQuantity
    requested_spacing: ScalarQuantity

    def __post_init__(self) -> None:
        values = (
            ("lower_bound", self.lower_bound),
            ("upper_bound", self.upper_bound),
            ("requested_spacing", self.requested_spacing),
        )
        for name, value in values:
            if not isinstance(value, ScalarQuantity):
                raise TypeError(f"{name} must be ScalarQuantity")
        if isinstance(self.lower_bound.unit, Unitless) != isinstance(
            self.upper_bound.unit, Unitless
        ) or isinstance(self.lower_bound.unit, Unitless) != isinstance(
            self.requested_spacing.unit, Unitless
        ):
            raise ValueError(
                "grid coordinates must not mix unitless and physical units"
            )
        for name, value in values[1:]:
            if not MODEL_SYSTEM_UNIT_CONVERTER.compatible(
                value.unit,
                self.lower_bound.unit,
            ):
                raise ValueError(f"{name} has incompatible coordinate dimensions")
        if self.upper_magnitude <= self.lower_magnitude:
            raise ValueError("upper_bound must be greater than lower_bound")
        if self.spacing_magnitude <= 0.0:
            raise ValueError("requested_spacing must be positive")
        interval_count = (
            self.upper_magnitude - self.lower_magnitude
        ) / self.spacing_magnitude
        rounded = int(round(interval_count))
        if not np.isclose(interval_count, rounded, rtol=0.0, atol=1.0e-12):
            raise ValueError("requested_spacing must divide the grid domain")
        if rounded < 1:
            raise ValueError("grid must contain at least one interval")

    @property
    def coordinate_unit(self) -> ModelSystemUnit:
        """Return the canonical unit selected by the lower bound."""
        return self.lower_bound.unit

    @property
    def lower_magnitude(self) -> float:
        """Return the lower coordinate in the canonical coordinate unit."""
        return self.lower_bound.magnitude

    @property
    def upper_magnitude(self) -> float:
        """Return the upper coordinate in the canonical coordinate unit."""
        return MODEL_SYSTEM_UNIT_CONVERTER.convert_scalar(
            self.upper_bound,
            self.coordinate_unit,
        ).magnitude

    @property
    def spacing_magnitude(self) -> float:
        """Return requested spacing in the canonical coordinate unit."""
        return MODEL_SYSTEM_UNIT_CONVERTER.convert_scalar(
            self.requested_spacing,
            self.coordinate_unit,
        ).magnitude

    @property
    def interval_count(self) -> int:
        """Return the exact number of uniform intervals."""
        return int(
            round(
                (self.upper_magnitude - self.lower_magnitude) / self.spacing_magnitude
            )
        )

    @property
    def point_count(self) -> int:
        """Return the number of grid points including both boundaries."""
        return self.interval_count + 1

    @property
    def interior_point_count(self) -> int:
        """Return the number of points excluding both boundaries."""
        return max(self.interval_count - 1, 0)

    @property
    def spacing(self) -> ScalarQuantity:
        """Return realized uniform spacing in the canonical unit."""
        magnitude = (self.upper_magnitude - self.lower_magnitude) / self.interval_count
        return ScalarQuantity(magnitude, self.coordinate_unit)

    def coordinates(self) -> VectorQuantity:
        """Return all ordered coordinates, including both boundaries."""
        values = self.lower_magnitude + self.spacing.magnitude * np.arange(
            self.point_count,
            dtype=np.float64,
        )
        return VectorQuantity(values, self.coordinate_unit)

    def interior_coordinates(self) -> VectorQuantity:
        """Return ordered coordinates excluding both boundary points."""
        return VectorQuantity(
            self.coordinates().magnitude[1:-1],
            self.coordinate_unit,
        )


@dataclass(frozen=True, slots=True)
class UniformCartesianGrid2D:
    """Represent a two-axis Cartesian tensor-product grid.

    The first and second axes must use compatible coordinate dimensions.
    Array representations use NumPy ``ij`` indexing and row-major ``C``
    flattening.
    """

    first_axis: UniformCartesianGrid1D
    second_axis: UniformCartesianGrid1D

    def __post_init__(self) -> None:
        if not isinstance(self.first_axis, UniformCartesianGrid1D):
            raise TypeError("first_axis must be UniformCartesianGrid1D")
        if not isinstance(self.second_axis, UniformCartesianGrid1D):
            raise TypeError("second_axis must be UniformCartesianGrid1D")
        if not MODEL_SYSTEM_UNIT_CONVERTER.compatible(
            self.first_axis.coordinate_unit,
            self.second_axis.coordinate_unit,
        ):
            raise ValueError("Cartesian axes must have compatible coordinate units")

    @property
    def shape(self) -> tuple[int, int]:
        """Return the tensor-product point shape including boundaries."""
        return (
            self.first_axis.point_count,
            self.second_axis.point_count,
        )

    @property
    def interior_shape(self) -> tuple[int, int]:
        """Return the tensor-product shape excluding every boundary face."""
        return (
            self.first_axis.interior_point_count,
            self.second_axis.interior_point_count,
        )

    @property
    def indexing(self) -> str:
        """Return the fixed NumPy mesh indexing convention ``ij``."""
        return "ij"

    @property
    def flattening_order(self) -> str:
        """Return the fixed row-major flattening order ``C``."""
        return "C"

    def coordinate_axes(self) -> tuple[VectorQuantity, VectorQuantity]:
        """Return the two ordered one-dimensional coordinate axes."""
        return (
            self.first_axis.coordinates(),
            self.second_axis.coordinates(),
        )
