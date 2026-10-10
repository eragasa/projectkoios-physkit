r"""Bloch-periodic Laplacians on nonorthogonal primitive-cell grids.

For Cartesian position :math:`\mathbf r=A\mathbf s`, the direct metric and its
inverse are

.. math::

    g=A^{\mathsf T}A,
    \qquad
    g^{-1}=\left(A^{\mathsf T}A\right)^{-1}.

Because the primitive basis is constant, the Cartesian Laplacian in fractional
coordinates is

.. math::

    \nabla_{\mathbf r}^{2}
    =\sum_{i,j}g^{ij}\frac{\partial^2}{\partial s_i\partial s_j}.

This module represents that geometric operator with centered second-order
finite differences and quotient-seam Bloch phases. It deliberately does not
apply kinetic-energy or effective-mass scaling.
"""

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
    r"""Request a centered geometric Laplacian in quotient-seam gauge.

    Parameters
    ----------
    grid
        Endpoint-excluded 2D or 3D fractional grid carrying direct-metric and
        basis-unit context.
    twist
        Boundary-twist lift :math:`\boldsymbol\kappa` measured in turns along
        the direct primitive axes. The component count must match the grid
        dimension. Components need not be reduced modulo one.

    Notes
    -----
    The requested Bloch convention is

    .. math::

        \psi(\mathbf s+\mathbf q)
        =\exp\!\left(2\pi i\,\boldsymbol\kappa\cdot\mathbf q\right)
        \psi(\mathbf s),
        \qquad \mathbf q\in\mathbb Z^d.

    The constructor retains the original lift while using its reduced
    representative for numerically bounded seam phases.

    Raises
    ------
    TypeError
        If ``grid`` or ``twist`` has the wrong nominal type.
    ValueError
        If the twist and grid dimensions differ.
    """

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
    r"""Correlate a represented Laplacian with its grid and twist context.

    Parameters
    ----------
    request
        Complete validated construction request.
    twist_fiber
        Reduction of the requested lift, explicitly tagged with
        ``TwistGaugeRepresentation.QUOTIENT_SEAM``.
    represented_laplacian
        Canonical immutable complex CSR matrix. Its square shape is the grid
        point count and its unit is unitless or inverse-length squared.

    Notes
    -----
    The matrix represents :math:`\nabla^2` rather than
    :math:`-\hbar^2\nabla^2/(2m)`. Consumers must compose any kinetic,
    constitutive, or potential terms in their owning model layer.

    Raises
    ------
    TypeError
        If a component has the wrong nominal type.
    ValueError
        If request, twist-fiber, matrix-shape, gauge, or unit metadata are
        inconsistent.
    """

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

        # Retain both the original lift and its exact reduced-seam context so a
        # caller never has to reconstruct which Bloch convention was used.
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

    For inverse metric entries :math:`g^{ij}` and fractional spacings
    :math:`h_i`, the represented operator is

    .. math::

        \Delta_h
        =\sum_i g^{ii}D_{ii}
        +2\sum_{i<j}g^{ij}D_{ij},

    with pure centered derivatives

    .. math::

        (D_{ii}f)_{\mathbf n}
        =\frac{f_{\mathbf n+\mathbf e_i}
        -2f_{\mathbf n}
        +f_{\mathbf n-\mathbf e_i}}{h_i^2},

    and mixed centered derivatives

    .. math::

        (D_{ij}f)_{\mathbf n}
        =\frac{
        f_{\mathbf n+\mathbf e_i+\mathbf e_j}
        -f_{\mathbf n+\mathbf e_i-\mathbf e_j}
        -f_{\mathbf n-\mathbf e_i+\mathbf e_j}
        +f_{\mathbf n-\mathbf e_i-\mathbf e_j}}
        {4h_i h_j}.

    If a raw target reduces according to

    .. math::

        \mathbf n+\boldsymbol\delta
        =\bar{\mathbf n}+N\mathbf q,

    its matrix contribution is multiplied by

    .. math::

        \exp\!\left(2\pi i\,\boldsymbol\kappa\cdot\mathbf q\right)
        =\prod_a\exp\!\left(2\pi i\kappa_a q_a\right).

    This product supplies the correct phase on simultaneous multi-axis wraps.
    """

    __slots__ = ()

    def action(
        self,
        *,
        request: BlochPeriodicLaplacianConstructionRequest,
    ) -> BlochPeriodicLaplacian:
        r"""Construct the canonical sparse geometric Laplacian.

        Parameters
        ----------
        request
            Validated grid and boundary-twist request.

        Returns
        -------
        BlochPeriodicLaplacian
            Correlated request, quotient-seam twist fiber, and canonical
            immutable complex CSR representation.

        Raises
        ------
        TypeError
            If ``request`` is not a
            `BlochPeriodicLaplacianConstructionRequest`.
        ValueError
            If a derived finite-difference coefficient is nonfinite. Component
            constructors can also reject inconsistent represented metadata.

        Notes
        -----
        Assembly remains sparse throughout. At most nine stencil locations per
        row occur in 2D and nineteen in 3D before zero mixed coefficients are
        omitted and duplicate COO entries are canonicalized.
        """
        if type(request) is not BlochPeriodicLaplacianConstructionRequest:
            raise TypeError(
                "request must be BlochPeriodicLaplacianConstructionRequest"
            )

        grid = request.grid
        dimension = grid.dimension.value
        extents = tuple(grid.domain.extents)
        spacings = grid.fractional_spacings
        inverse_metric = grid.metric.inverse_metric_tensor

        # Quotient-seam phases depend only on kappa modulo integer shifts. Keep
        # the full reduction in the result, but evaluate bounded representatives.
        reduction = BoundaryTwistReducer().execute(request.twist)
        representative = reduction.representative.turns
        twist_fiber = TwistFiber(
            reduction,
            TwistGaugeRepresentation.QUOTIENT_SEAM,
        )

        # Metric coefficients are position-independent, so build the complete
        # pure-plus-mixed stencil once and reuse it for every grid row.
        stencil = self._stencil(
            inverse_metric=inverse_metric,
            spacings=spacings,
            dimension=dimension,
        )

        rows: list[int] = []
        columns: list[int] = []
        values: list[complex] = []
        coordinate_ranges = tuple(range(extent) for extent in extents)

        # itertools.product varies its final iterable fastest, matching the
        # FinitePeriodicDomain C-order contract and making enumerate() the row.
        for row, source in enumerate(product(*coordinate_ranges)):
            for displacement, coefficient in stencil:
                wrapped: list[int] = []
                quotients: list[int] = []
                for axis in range(dimension):
                    # Euclidean divmod returns both the in-domain coordinate and
                    # the signed number of primitive-cell boundaries crossed.
                    quotient, coordinate = divmod(
                        source[axis] + displacement[axis],
                        extents[axis],
                    )
                    quotients.append(quotient)
                    wrapped.append(coordinate)

                # Summing all quotient components before exponentiation is
                # equivalent to multiplying one seam phase per wrapped axis.
                phase_turns = math.fsum(
                    representative[axis] * quotients[axis]
                    for axis in range(dimension)
                )
                phase = cmath.exp(2.0j * math.pi * phase_turns)
                rows.append(row)
                columns.append(self._c_order_index(wrapped, extents))
                values.append(coefficient * phase)

        point_count = grid.point_count
        # COO naturally accumulates stencil collisions. The quantity boundary
        # converts it to sorted, duplicate-free, explicit-zero-free CSR storage.
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
        r"""Return the constant pure- and mixed-derivative stencil.

        Parameters
        ----------
        inverse_metric
            Finite symmetric matrix containing :math:`g^{ij}`.
        spacings
            Fractional spacings :math:`h_i=1/N_i`.
        dimension
            Spatial dimension, already validated as two or three.

        Returns
        -------
        tuple
            Ordered ``(displacement, coefficient)`` pairs. Each displacement
            is an integer fractional-grid offset.

        Raises
        ------
        ValueError
            If any assembled coefficient is nonfinite.

        Notes
        -----
        The center coefficient is

        .. math::

            -2\sum_i\frac{g^{ii}}{h_i^2}.

        Axial neighbors carry :math:`g^{ii}/h_i^2`. For ``i < j``, a mixed
        neighbor with signs :math:`\sigma_i,\sigma_j\in\{-1,+1\}` carries

        .. math::

            \frac{g^{ij}\sigma_i\sigma_j}{2h_i h_j},

        where the denominator includes the factor two arising from
        :math:`2g^{ij}D_{ij}`.
        """
        zero = (0,) * dimension
        center = math.fsum(
            -2.0 * float(inverse_metric[axis, axis]) / spacings[axis] ** 2
            for axis in range(dimension)
        )
        terms: list[tuple[tuple[int, ...], float]] = [(zero, center)]

        # Pure second derivatives contribute one equal real coefficient to the
        # positive and negative neighbor on each fractional axis.
        for axis in range(dimension):
            coefficient = (
                float(inverse_metric[axis, axis]) / spacings[axis] ** 2
            )
            for sign in (-1, 1):
                displacement = [0] * dimension
                displacement[axis] = sign
                terms.append((tuple(displacement), coefficient))

        # Off-diagonal inverse-metric entries require the four-corner mixed
        # stencil. Same-sign corners are positive; opposite-sign corners negate.
        for first_axis in range(dimension):
            for second_axis in range(first_axis + 1, dimension):
                scale = float(inverse_metric[first_axis, second_axis]) / (
                    2.0 * spacings[first_axis] * spacings[second_axis]
                )
                # Orthogonal axis pairs have no mixed term and need no zeros in
                # the sparse COO inventory.
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
        r"""Map one wrapped coordinate to its C-order linear index.

        Parameters
        ----------
        coordinates
            In-domain integer coordinates in fractional-axis order.
        extents
            Matching axis point counts.

        Returns
        -------
        int
            Last-axis-fastest index obtained recursively from
            :math:`I_{a+1}=I_aN_{a+1}+n_{a+1}`.
        """
        result = 0
        for coordinate, extent in zip(coordinates, extents, strict=True):
            result = result * extent + coordinate
        return result
