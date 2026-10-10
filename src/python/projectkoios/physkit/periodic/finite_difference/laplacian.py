r"""Bloch-periodic Laplacians on nonorthogonal primitive-cell grids."""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy import sparse

from projectkoios.physkit.core.actions import DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.finite_difference.grid import (
    PeriodicFiniteDifferenceGrid,
)
from projectkoios.physkit.periodic.lattice.boundary_phases import (
    BoundaryTwistLift,
    BoundaryTwistReducer,
    TwistFiber,
    TwistGaugeRepresentation,
)
from projectkoios.physkit.units import ComplexSparseMatrixQuantity


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class BlochPeriodicLaplacianConstructionRequest(DataObject):
    """Request a centered geometric Laplacian in quotient-seam gauge."""

    grid: PeriodicFiniteDifferenceGrid
    twist: BoundaryTwistLift

    def __post_init__(self) -> None:
        """Validate the grid and dimension-compatible twist lift."""
        if type(self.grid) is not PeriodicFiniteDifferenceGrid:
            raise TypeError("grid must be PeriodicFiniteDifferenceGrid")
        if type(self.twist) is not BoundaryTwistLift:
            raise TypeError("twist must be BoundaryTwistLift")
        if self.twist.dimension is not self.grid.dimension:
            raise ValueError("twist and grid dimensions must agree")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class BlochPeriodicLaplacian(ResultsObject):
    """Return a represented geometric Laplacian and its complete context."""

    request: BlochPeriodicLaplacianConstructionRequest
    twist_fiber: TwistFiber
    represented_laplacian: ComplexSparseMatrixQuantity

    def __post_init__(self) -> None:
        """Validate request, quotient-seam metadata, matrix shape, and unit."""
        if type(self.request) is not BlochPeriodicLaplacianConstructionRequest:
            raise TypeError(
                "request must be BlochPeriodicLaplacianConstructionRequest"
            )
        if type(self.twist_fiber) is not TwistFiber:
            raise TypeError("twist_fiber must be TwistFiber")
        if self.twist_fiber.lift != self.request.twist:
            raise ValueError("twist_fiber must retain the requested twist lift")
        if self.twist_fiber.dimension is not self.request.grid.dimension:
            raise ValueError("twist_fiber and grid dimensions must agree")
        if self.twist_fiber.gauge is not TwistGaugeRepresentation.QUOTIENT_SEAM:
            raise ValueError("twist_fiber must use quotient-seam gauge")
        if type(self.represented_laplacian) is not ComplexSparseMatrixQuantity:
            raise TypeError(
                "represented_laplacian must be ComplexSparseMatrixQuantity"
            )
        expected_shape = (
            self.request.grid.point_count,
            self.request.grid.point_count,
        )
        if self.represented_laplacian.shape != expected_shape:
            raise ValueError("represented_laplacian shape must match the grid")
        if self.represented_laplacian.unit != self.request.grid.laplacian_unit:
            raise ValueError("represented_laplacian unit must match the grid")


class BlochPeriodicLaplacianConstructor(
    DataObjectActionizer[
        BlochPeriodicLaplacianConstructionRequest,
        BlochPeriodicLaplacian,
    ]
):
    r"""Construct the metric Laplacian with Bloch-twisted periodic seams.

    For direct basis ``A``, fractional coordinates ``r=A s``, inverse metric
    ``g^{-1}``, and fractional spacings ``h_i``, the represented operator is

    .. math::

        \Delta_h = \sum_i g^{ii}D_{ii}
        + 2\sum_{i<j}g^{ij}D_{ij}.

    Every raw stencil target is reduced into the finite domain. A boundary
    quotient ``q`` contributes the quotient-seam factor
    ``exp(2 pi i kappa dot q)`` using the reduced twist representative.
    """

    __slots__ = ()

    def action(
        self,
        *,
        request: BlochPeriodicLaplacianConstructionRequest,
    ) -> BlochPeriodicLaplacian:
        """Return a canonical complex CSR representation of the Laplacian."""
        if type(request) is not BlochPeriodicLaplacianConstructionRequest:
            raise TypeError(
                "request must be BlochPeriodicLaplacianConstructionRequest"
            )

        grid = request.grid
        dimension = grid.dimension.value
        extents = tuple(grid.domain.extents)
        spacings = grid.fractional_spacings
        inverse_metric = grid.metric.inverse_metric_tensor
        reduction = BoundaryTwistReducer().execute(request.twist)
        representative = reduction.representative.turns
        twist_fiber = TwistFiber(
            reduction,
            TwistGaugeRepresentation.QUOTIENT_SEAM,
        )
        stencil = self._stencil(
            inverse_metric=inverse_metric,
            spacings=spacings,
            dimension=dimension,
        )

        rows: list[int] = []
        columns: list[int] = []
        values: list[complex] = []
        coordinate_ranges = tuple(range(extent) for extent in extents)
        for row, source in enumerate(product(*coordinate_ranges)):
            for displacement, coefficient in stencil:
                wrapped: list[int] = []
                quotients: list[int] = []
                for axis in range(dimension):
                    quotient, coordinate = divmod(
                        source[axis] + displacement[axis],
                        extents[axis],
                    )
                    quotients.append(quotient)
                    wrapped.append(coordinate)
                phase_turns = math.fsum(
                    representative[axis] * quotients[axis]
                    for axis in range(dimension)
                )
                phase = cmath.exp(2.0j * math.pi * phase_turns)
                rows.append(row)
                columns.append(self._c_order_index(wrapped, extents))
                values.append(coefficient * phase)

        point_count = grid.point_count
        assembled = sparse.coo_array(
            (
                np.asarray(values, dtype=np.complex128),
                (
                    np.asarray(rows, dtype=np.int64),
                    np.asarray(columns, dtype=np.int64),
                ),
            ),
            shape=(point_count, point_count),
            dtype=np.complex128,
        )
        return BlochPeriodicLaplacian(
            request=request,
            twist_fiber=twist_fiber,
            represented_laplacian=ComplexSparseMatrixQuantity.from_csr(
                assembled,
                grid.laplacian_unit,
            ),
        )

    @staticmethod
    def _stencil(
        *,
        inverse_metric: np.ndarray,
        spacings: tuple[float, ...],
        dimension: int,
    ) -> tuple[tuple[tuple[int, ...], float], ...]:
        """Return centered pure- and mixed-derivative stencil terms."""
        zero = (0,) * dimension
        center = math.fsum(
            -2.0 * float(inverse_metric[axis, axis]) / spacings[axis] ** 2
            for axis in range(dimension)
        )
        terms: list[tuple[tuple[int, ...], float]] = [(zero, center)]

        for axis in range(dimension):
            coefficient = (
                float(inverse_metric[axis, axis]) / spacings[axis] ** 2
            )
            for sign in (-1, 1):
                displacement = [0] * dimension
                displacement[axis] = sign
                terms.append((tuple(displacement), coefficient))

        for first_axis in range(dimension):
            for second_axis in range(first_axis + 1, dimension):
                scale = float(inverse_metric[first_axis, second_axis]) / (
                    2.0 * spacings[first_axis] * spacings[second_axis]
                )
                if scale == 0.0:
                    continue
                for first_sign in (-1, 1):
                    for second_sign in (-1, 1):
                        displacement = [0] * dimension
                        displacement[first_axis] = first_sign
                        displacement[second_axis] = second_sign
                        coefficient = scale * first_sign * second_sign
                        terms.append((tuple(displacement), coefficient))

        if any(not math.isfinite(coefficient) for _, coefficient in terms):
            raise ValueError("finite-difference stencil coefficients must be finite")
        return tuple(terms)

    @staticmethod
    def _c_order_index(coordinates: list[int], extents: tuple[int, ...]) -> int:
        """Return the last-axis-fastest linear index."""
        result = 0
        for coordinate, extent in zip(coordinates, extents, strict=True):
            result = result * extent + coordinate
        return result
