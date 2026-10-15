"""One-dimensional uniform-grid representations.

This module retains the legacy :class:`Grid1D` active-set representation and
adds :class:`UniformGrid1D`, an immutable geometry-only representation of a
closed interval.  ``UniformGrid1D`` includes both finite endpoints, uses
binary64 coordinates, and requires at least three points.  It deliberately
has no boundary-condition, active-index, restriction, or embedding behavior;
those homogeneous-Dirichlet semantics belong to
:mod:`projectkoios.physkit.discretization.state_space_1d`.

The grid spacing is ``(b - a) / (num_points - 1)``.  Coordinate arrays are
owned, C-contiguous copies, so callers cannot mutate retained state.  Invalid
scalar types raise :class:`TypeError`; finite, ordering, and size invariant
violations raise :class:`ValueError`.  This module defines geometry and makes
no physical-model, accuracy, convergence, validation, or uncertainty claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from numbers import Integral, Real

import numpy as np


class ActiveSetType1D(Enum):
    """Active index sets on a full closed reference grid."""

    ALL = auto()
    INTERIOR = auto()
    LEFT_CLOSED = auto()
    RIGHT_CLOSED = auto()
    LEFT_BOUNDARY = auto()
    RIGHT_BOUNDARY = auto()
    BOUNDARY = auto()


@dataclass(frozen=True)
class Grid1D:
    """Uniform 1D reference grid on the closed interval [a,b].

    The full coordinate grid always includes both endpoints:

        x = np.linspace(a, b, N)

    The active set selects which grid points are used as degrees of freedom.
    """

    a: float
    b: float
    N: int
    active_type: ActiveSetType1D = ActiveSetType1D.ALL

    def __post_init__(self) -> None:
        if self.b <= self.a:
            raise ValueError("Require b > a.")

        if self.N < 2:
            raise ValueError("Require N >= 2.")

    @property
    def L(self) -> float:
        return self.b - self.a

    @property
    def dx(self) -> float:
        return self.L / (self.N - 1)

    @property
    def x(self) -> np.ndarray:
        return np.linspace(self.a, self.b, self.N)

    @property
    def active_indices(self) -> np.ndarray:
        if self.active_type is ActiveSetType1D.ALL:
            return np.arange(0, self.N)

        if self.active_type is ActiveSetType1D.INTERIOR:
            return np.arange(1, self.N - 1)

        if self.active_type is ActiveSetType1D.LEFT_CLOSED:
            return np.arange(0, self.N - 1)

        if self.active_type is ActiveSetType1D.RIGHT_CLOSED:
            return np.arange(1, self.N)

        if self.active_type is ActiveSetType1D.LEFT_BOUNDARY:
            return np.array([0])

        if self.active_type is ActiveSetType1D.RIGHT_BOUNDARY:
            return np.array([self.N - 1])

        if self.active_type is ActiveSetType1D.BOUNDARY:
            return np.array([0, self.N - 1])

        raise ValueError(f"Unknown active_type: {self.active_type}")

    @property
    def x_active(self) -> np.ndarray:
        return self.x[self.active_indices]


class UniformGrid1D:
    """Immutable uniform geometry on a closed one-dimensional interval.

    Parameters
    ----------
    a : numbers.Real
        Finite left endpoint.
    b : numbers.Real
        Finite right endpoint, strictly greater than ``a``.
    num_points : numbers.Integral
        Number of endpoint-inclusive points; at least three.

    Attributes
    ----------
    a : float
        Left endpoint stored as a built-in ``float``.
    b : float
        Right endpoint stored as a built-in ``float``.
    num_points : int
        Number of grid points stored as a built-in ``int``.
    length : float
        Interval length ``b - a``.
    spacing : float
        Uniform spacing ``(b - a) / (num_points - 1)``.
    coordinates : numpy.ndarray
        Owned C-contiguous ``float64`` endpoint-inclusive coordinates.

    Raises
    ------
    TypeError
        If an endpoint is Boolean or not a Python/NumPy real scalar, or if
        ``num_points`` is Boolean or not a Python/NumPy integral scalar.
    ValueError
        If an endpoint is nonfinite, ``b <= a``, or ``num_points < 3``.

    Notes
    -----
    This object owns geometry only.  Homogeneous-Dirichlet boundaries and
    active-vector interpretation are provided by
    :class:`~projectkoios.physkit.discretization.HomogeneousDirichletStateSpace1D`.
    Returned coordinates never alias retained storage.

    Examples
    --------
    >>> grid = UniformGrid1D(0.0, 1.0, 5)
    >>> grid.spacing
    0.25
    >>> grid.coordinates
    array([0.  , 0.25, 0.5 , 0.75, 1.  ])

    See Also
    --------
    projectkoios.physkit.discretization.HomogeneousDirichletStateSpace1D
        State-space semantics for this geometry.
    """

    __slots__ = ("_a", "_b", "_num_points", "_coordinates", "_frozen")

    def __init__(self, a: Real, b: Real, num_points: Integral) -> None:
        if isinstance(a, (bool, np.bool_)) or not isinstance(a, Real):
            raise TypeError("a must be a non-Boolean real scalar.")
        if isinstance(b, (bool, np.bool_)) or not isinstance(b, Real):
            raise TypeError("b must be a non-Boolean real scalar.")
        if isinstance(num_points, (bool, np.bool_)) or not isinstance(
            num_points, Integral
        ):
            raise TypeError("num_points must be a non-Boolean integral scalar.")

        try:
            a_value = float(a)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("a cannot be canonicalized as a finite float.") from exc
        try:
            b_value = float(b)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("b cannot be canonicalized as a finite float.") from exc
        point_count = int(num_points)
        if not np.isfinite(a_value):
            raise ValueError("a must be finite.")
        if not np.isfinite(b_value):
            raise ValueError("b must be finite.")
        if b_value <= a_value:
            raise ValueError("b must be greater than a.")
        if point_count < 3:
            raise ValueError("num_points must be at least 3.")

        object.__setattr__(self, "_a", a_value)
        object.__setattr__(self, "_b", b_value)
        object.__setattr__(self, "_num_points", point_count)
        coordinates = np.linspace(a_value, b_value, point_count, dtype=np.float64)
        object.__setattr__(self, "_coordinates", coordinates)
        object.__setattr__(self, "_frozen", True)

    def __setattr__(self, name: str, value: object) -> None:
        """Reject reassignment after construction.

        Raises
        ------
        AttributeError
            Always, after the instance has been successfully constructed.
        """
        if getattr(self, "_frozen", False):
            raise AttributeError(f"{type(self).__name__} is immutable.")
        object.__setattr__(self, name, value)

    @property
    def a(self) -> float:
        """Return the left endpoint.

        Returns
        -------
        float
            Finite left endpoint.
        """
        return self._a

    @property
    def b(self) -> float:
        """Return the right endpoint.

        Returns
        -------
        float
            Finite right endpoint.
        """
        return self._b

    @property
    def num_points(self) -> int:
        """Return the number of endpoint-inclusive points.

        Returns
        -------
        int
            Grid size, at least three.
        """
        return self._num_points

    @property
    def length(self) -> float:
        """Return the interval length.

        Returns
        -------
        float
            ``b - a``.
        """
        return self._b - self._a

    @property
    def spacing(self) -> float:
        """Return the uniform closed-grid spacing.

        Returns
        -------
        float
            ``(b - a) / (num_points - 1)``.
        """
        return self.length / (self._num_points - 1)

    @property
    def coordinates(self) -> np.ndarray:
        """Return an owned copy of the full coordinates.

        Returns
        -------
        numpy.ndarray
            New C-contiguous ``float64`` array of shape ``(num_points,)``.
        """
        return self._coordinates.copy(order="C")
