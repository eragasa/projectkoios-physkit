"""Homogeneous-Dirichlet state-space semantics on a uniform 1D grid.

The module separates immutable geometry from representation semantics.  A
:class:`HomogeneousDirichletStateSpace1D` contains a closed
:class:`~projectkoios.physkit.discretization.grid_1d.UniformGrid1D`, prescribes
exact zero endpoint values, and interprets the ``N - 2`` interior samples as a real or
complex active vector.  Restriction discards endpoints without checking their
values; embedding supplies exact zeros.

Accepted vectors are finite, numeric, one-dimensional arrays of the exact
required shape.  Boolean, object, string, and ragged inputs are rejected.
Real results are owned C-contiguous ``float64`` arrays and complex results are
owned C-contiguous ``complex128`` arrays.  Wrong semantic types raise
:class:`TypeError`; rank, shape, and finiteness violations raise
:class:`ValueError`.  These representation rules make no physical-model,
convergence, validation, or uncertainty claim.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from .grid_1d import UniformGrid1D


def _validated_numeric_vector(
    values: ArrayLike, expected_length: int, name: str
) -> np.ndarray:
    """Canonicalize one finite real or complex vector without retaining it."""
    try:
        uncoerced = np.asarray(values, dtype=object)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a non-ragged numeric vector.") from exc
    if any(isinstance(value, (bool, np.bool_)) for value in uncoerced.flat):
        raise TypeError(f"{name} must not contain Boolean values.")

    try:
        array = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a non-ragged numeric vector.") from exc

    if array.dtype.kind not in "iufc":
        raise TypeError(f"{name} must contain non-Boolean real or complex numbers.")
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional.")
    if array.shape != (expected_length,):
        raise ValueError(f"{name} must have shape ({expected_length},).")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values.")

    dtype = np.complex128 if np.iscomplexobj(array) else np.float64
    result = np.array(array, dtype=dtype, order="C", copy=True)
    if not np.all(np.isfinite(result)):
        raise ValueError(
            f"{name} cannot be represented by finite {dtype.__name__} values."
        )
    return result


class HomogeneousDirichletStateSpace1D:
    """Immutable active-vector semantics for zero endpoint values.

    Parameters
    ----------
    grid : UniformGrid1D
        Immutable closed, endpoint-inclusive geometry.

    Attributes
    ----------
    grid : UniformGrid1D
        The exact contained grid instance.
    boundary_values : tuple[float, float]
        Exact prescribed endpoint values ``(0.0, 0.0)``.
    active_indices : numpy.ndarray
        Owned platform-integer indices ``1`` through ``N - 2``.
    active_coordinates : numpy.ndarray
        Owned C-contiguous ``float64`` interior coordinates.
    dimension : int
        Active-vector dimension ``N - 2``.
    semantic_identity : tuple
        Exact compatibility identity including geometry, endpoints, active
        convention, boundary convention, and real/complex interpretation.

    Raises
    ------
    TypeError
        If ``grid`` is not a :class:`UniformGrid1D`.

    Notes
    -----
    Compatibility is equality of :attr:`semantic_identity`, not object
    identity, shape alone, or approximate coordinate equality.  Arrays
    returned by this object do not alias retained storage.

    Examples
    --------
    >>> from projectkoios.physkit.discretization import UniformGrid1D
    >>> space = HomogeneousDirichletStateSpace1D(UniformGrid1D(0.0, 1.0, 5))
    >>> space.dimension
    3
    >>> space.embed([1.0, 2.0, 3.0])
    array([0., 1., 2., 3., 0.])

    See Also
    --------
    UniformGrid1D
        The contained geometry-only grid.
    projectkoios.physkit.operators.FiniteDifferenceLaplacian1D
        An operator whose domain and codomain are this state space.
    """

    __slots__ = (
        "_grid",
        "_active_indices",
        "_active_coordinates",
        "_semantic_identity",
        "_frozen",
    )

    def __init__(self, grid: UniformGrid1D) -> None:
        if not isinstance(grid, UniformGrid1D):
            raise TypeError("grid must be a UniformGrid1D.")

        active_indices = np.arange(1, grid.num_points - 1, dtype=np.intp)
        active_coordinates = np.array(
            grid.coordinates[active_indices], dtype=np.float64, order="C", copy=True
        )
        identity = (
            "HomogeneousDirichletStateSpace1D",
            (
                "UniformGrid1D",
                float(grid.a),
                float(grid.b),
                int(grid.num_points),
                "closed-endpoint-inclusive",
            ),
            "active-interior-indices-1-through-N-minus-2",
            ("homogeneous-dirichlet", 0.0, 0.0),
            "real-or-complex-active-vectors",
        )

        object.__setattr__(self, "_grid", grid)
        object.__setattr__(self, "_active_indices", active_indices)
        object.__setattr__(self, "_active_coordinates", active_coordinates)
        object.__setattr__(self, "_semantic_identity", identity)
        object.__setattr__(self, "_frozen", True)

    def __setattr__(self, name: str, value: Any) -> None:
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
    def grid(self) -> UniformGrid1D:
        """Return the contained immutable grid.

        Returns
        -------
        UniformGrid1D
            Exact grid instance supplied at construction.
        """
        return self._grid

    @property
    def boundary_values(self) -> tuple[float, float]:
        """Return the prescribed endpoint values.

        Returns
        -------
        tuple[float, float]
            ``(0.0, 0.0)``.
        """
        return (0.0, 0.0)

    @property
    def active_indices(self) -> np.ndarray:
        """Return owned interior grid indices.

        Returns
        -------
        numpy.ndarray
            New platform-integer array of shape ``(dimension,)``.
        """
        return self._active_indices.copy(order="C")

    @property
    def active_coordinates(self) -> np.ndarray:
        """Return owned interior coordinates.

        Returns
        -------
        numpy.ndarray
            New C-contiguous ``float64`` array of shape ``(dimension,)``.
        """
        return self._active_coordinates.copy(order="C")

    @property
    def dimension(self) -> int:
        """Return the active-vector dimension.

        Returns
        -------
        int
            ``grid.num_points - 2``.
        """
        return self._grid.num_points - 2

    @property
    def semantic_identity(self) -> tuple:
        """Return the exact representation-compatibility identity.

        Returns
        -------
        tuple
            Immutable value tuple specified by the capability contract.

        Notes
        -----
        Equality denotes representation compatibility only; it does not
        establish cross-grid physical equivalence.
        """
        return self._semantic_identity

    def restrict(self, full_state: ArrayLike) -> np.ndarray:
        """Restrict a full sampled vector to its interior entries.

        Parameters
        ----------
        full_state : numpy.typing.ArrayLike
            Finite numeric one-dimensional vector of shape
            ``(grid.num_points,)``.

        Returns
        -------
        numpy.ndarray
            New owned C-contiguous active vector.  Real input produces
            ``float64`` and complex input produces ``complex128``.

        Raises
        ------
        TypeError
            If the input is Boolean, object, string, ragged, or otherwise
            nonnumeric.
        ValueError
            If rank, shape, finiteness, or canonical representation is
            invalid.

        Notes
        -----
        Endpoint values are discarded without validation or correction.

        See Also
        --------
        embed
            Insert an active vector between exact zero endpoints.
        """
        full = _validated_numeric_vector(
            full_state, self._grid.num_points, "full_state"
        )
        return np.array(full[1:-1], dtype=full.dtype, order="C", copy=True)

    def embed(self, active_state: ArrayLike) -> np.ndarray:
        """Embed an active vector between exact zero endpoints.

        Parameters
        ----------
        active_state : numpy.typing.ArrayLike
            Finite numeric one-dimensional vector of shape ``(dimension,)``.

        Returns
        -------
        numpy.ndarray
            New owned C-contiguous full vector.  Real input produces
            ``float64`` and complex input produces ``complex128``.

        Raises
        ------
        TypeError
            If the input is Boolean, object, string, ragged, or otherwise
            nonnumeric.
        ValueError
            If rank, shape, finiteness, or canonical representation is
            invalid.

        Notes
        -----
        The inserted endpoint values are exact zeros of the result dtype.

        See Also
        --------
        restrict
            Project a full sampled vector onto its interior entries.
        """
        active = _validated_numeric_vector(active_state, self.dimension, "active_state")
        full = np.zeros(self._grid.num_points, dtype=active.dtype, order="C")
        full[1:-1] = active
        return full
