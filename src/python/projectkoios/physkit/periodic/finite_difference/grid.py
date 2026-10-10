"""Endpoint-excluded fractional grids on two- and three-dimensional cells."""

from __future__ import annotations

from dataclasses import dataclass

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.periodic.lattice.finite_domain import (
    FinitePeriodicDomain,
    LatticeDimension,
    LatticeIntegerComponents,
)
from projectkoios.physkit.periodic.lattice.metric import DirectLatticeMetric
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    PhysicalUnit,
    Unitless,
)

_METRE = PhysicalUnit("meter")

type FractionalGridSpacings = tuple[float, float] | tuple[float, float, float]
"""Closed fractional-spacing representation for two or three dimensions."""


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicFiniteDifferenceGrid(DataObject):
    r"""Define an endpoint-excluded fractional grid on one primitive cell.

    The direct-basis magnitudes owned by ``metric`` are interpreted in
    ``direct_basis_unit``. Grid sites have fractional coordinates
    ``s_i = n_i / N_i`` for ``n_i`` from zero through ``N_i - 1``. The finite
    domain fixes C-order flattening, with the last fractional axis varying
    fastest.
    """

    metric: DirectLatticeMetric
    domain: FinitePeriodicDomain
    direct_basis_unit: PhysicalUnit | Unitless

    def __post_init__(self) -> None:
        """Validate geometry, finite extents, and the basis-unit declaration."""
        if type(self.metric) is not DirectLatticeMetric:
            raise TypeError("metric must be DirectLatticeMetric")
        if type(self.domain) is not FinitePeriodicDomain:
            raise TypeError("domain must be FinitePeriodicDomain")
        if self.domain.dimension not in (
            LatticeDimension.TWO,
            LatticeDimension.THREE,
        ):
            raise ValueError("domain dimension must be two or three")
        if self.metric.dimension != self.domain.dimension.value:
            raise ValueError("metric and domain dimensions must agree")
        if any(extent < 3 for extent in self.domain.extents):
            raise ValueError("every periodic grid extent must be at least three")
        if not isinstance(self.direct_basis_unit, PhysicalUnit | Unitless):
            raise TypeError("direct_basis_unit must be PhysicalUnit or Unitless")
        if isinstance(self.direct_basis_unit, PhysicalUnit) and not (
            MODEL_SYSTEM_UNIT_CONVERTER.compatible(
                self.direct_basis_unit,
                _METRE,
            )
        ):
            raise ValueError("direct_basis_unit must be compatible with length")

    @property
    def dimension(self) -> LatticeDimension:
        """Return the exact supported spatial dimension."""
        return self.domain.dimension

    @property
    def shape(self) -> LatticeIntegerComponents:
        """Return endpoint-excluded point counts in fractional-axis order."""
        return self.domain.extents

    @property
    def point_count(self) -> int:
        """Return the tensor-product grid-point count."""
        return self.domain.cell_count

    @property
    def fractional_spacings(self) -> FractionalGridSpacings:
        """Return dimensionless spacings ``h_i = 1 / N_i``."""
        values = tuple(1.0 / extent for extent in self.domain.extents)
        if self.dimension is LatticeDimension.TWO:
            return (values[0], values[1])
        return (values[0], values[1], values[2])

    @property
    def laplacian_unit(self) -> PhysicalUnit | Unitless:
        """Return the represented geometric Laplacian unit."""
        if isinstance(self.direct_basis_unit, Unitless):
            return Unitless()
        return PhysicalUnit(
            f"1 / ({self.direct_basis_unit.expression}) ** 2"
        )
