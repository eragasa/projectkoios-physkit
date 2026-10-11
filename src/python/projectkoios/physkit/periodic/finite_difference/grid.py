r"""Endpoint-excluded fractional grids on periodic primitive cells.

A primitive-cell point is parameterized by fractional coordinates
:math:`\mathbf s` through

.. math::

    \mathbf r = A\mathbf s,
    \qquad
    s_i = \frac{n_i}{N_i},
    \qquad
    n_i\in\{0,\ldots,N_i-1\}.

The repeated endpoint :math:`s_i=1` is excluded because it is identified with
:math:`s_i=0` by the periodic boundary condition. This module owns grid
metadata only; represented differential operators live in :mod:`laplacian`.
"""

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
    r"""Define an endpoint-excluded grid on one periodic primitive cell.

    Parameters
    ----------
    metric
        Metric geometry derived from an intrinsic ``2 x 2`` or full ``3 x 3``
        direct primitive basis. The basis vectors are the columns of
        :math:`A`.
    domain
        Finite periodic index domain. Its extents :math:`(N_0,\ldots,N_{d-1})`
        are the point counts along the fractional axes. Every extent must be at
        least three so that the positive and negative centered-stencil
        neighbors remain distinct.
    direct_basis_unit
        Unit in which the numerical entries of ``metric.primitive_basis`` are
        represented. A physical value must be length-compatible; ``Unitless``
        declares a nondimensional basis.

    Notes
    -----
    Grid points satisfy

    .. math::

        \mathbf s_{\mathbf n}
        = \left(\frac{n_0}{N_0},\ldots,
        \frac{n_{d-1}}{N_{d-1}}\right),
        \qquad
        h_i=\frac{1}{N_i}.

    The endpoint is omitted on every axis. Array states use NumPy C-order
    flattening, so the final fractional axis varies fastest. The unit field is
    declarative: this record does not rescale basis magnitudes.

    Raises
    ------
    TypeError
        If a composed object or unit has the wrong nominal type.
    ValueError
        If the dimension is not two or three, component dimensions disagree,
        an extent is smaller than three, or a physical unit is not a length.
    """

    metric: DirectLatticeMetric
    domain: FinitePeriodicDomain
    direct_basis_unit: PhysicalUnit | Unitless

    def __post_init__(self) -> None:
        """Validate geometry, finite extents, and the basis-unit declaration."""
        # Require the nominal metric and domain owners rather than accepting
        # structurally similar objects with unknown conventions.
        if type(self.metric) is not DirectLatticeMetric:
            raise TypeError("metric must be DirectLatticeMetric")
        if type(self.domain) is not FinitePeriodicDomain:
            raise TypeError("domain must be FinitePeriodicDomain")

        # The shared stencil is deliberately closed to the agreed 2D/3D scope.
        if self.domain.dimension not in (
            LatticeDimension.TWO,
            LatticeDimension.THREE,
        ):
            raise ValueError("domain dimension must be two or three")
        if self.metric.dimension != self.domain.dimension.value:
            raise ValueError("metric and domain dimensions must agree")

        # With N_i >= 3, offsets +e_i and -e_i address distinct grid sites.
        # This prevents centered pure and mixed derivatives from degenerating.
        if any(extent < 3 for extent in self.domain.extents):
            raise ValueError("every periodic grid extent must be at least three")

        # Units describe the numerical basis scale; physical bases must carry a
        # true length unit so the Laplacian can expose inverse-length squared.
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
        """Return the exact supported spatial dimension.

        Returns
        -------
        LatticeDimension
            ``LatticeDimension.TWO`` or ``LatticeDimension.THREE``.
        """
        return self.domain.dimension

    @property
    def shape(self) -> LatticeIntegerComponents:
        """Return point counts in fractional-axis order.

        Returns
        -------
        LatticeIntegerComponents
            Two or three endpoint-excluded axis extents. The final extent is
            the fastest-varying axis under the fixed C-order convention.
        """
        return self.domain.extents

    @property
    def point_count(self) -> int:
        r"""Return the tensor-product grid-point count.

        Returns
        -------
        int
            The exact product :math:`\prod_i N_i`.
        """
        return self.domain.cell_count

    @property
    def fractional_spacings(self) -> FractionalGridSpacings:
        r"""Return dimensionless fractional-coordinate spacings.

        Returns
        -------
        FractionalGridSpacings
            Ordered values :math:`h_i=1/N_i`, one for each fractional axis.
        """
        values = tuple(1.0 / extent for extent in self.domain.extents)
        # Close the return type explicitly instead of leaking an arbitrary tuple.
        if self.dimension is LatticeDimension.TWO:
            return (values[0], values[1])
        return (values[0], values[1], values[2])

    @property
    def laplacian_unit(self) -> PhysicalUnit | Unitless:
        r"""Return the represented geometric Laplacian unit.

        Returns
        -------
        PhysicalUnit or Unitless
            ``Unitless`` for a nondimensional basis; otherwise
            :math:`[A]^{-2}`, the inverse square of the declared basis unit.
        """
        if isinstance(self.direct_basis_unit, Unitless):
            return Unitless()
        return PhysicalUnit(f"1 / ({self.direct_basis_unit.expression}) ** 2")
