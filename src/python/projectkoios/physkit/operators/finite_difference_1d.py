"""Centered finite-difference second derivative on a uniform 1D grid.

This module implements :class:`FiniteDifferenceLaplacian1D`, the eager CSR
representation of the positive second-derivative convention ``+d**2/dx**2``
on a homogeneous-Dirichlet active state space.  For spacing ``h`` and
``M = N - 2`` active values, the matrix is
``h**-2 * tridiag(1, -2, 1)`` with shape ``(M, M)``.  It is real symmetric
negative definite; for ``N = 3`` it is exactly the one-by-one matrix
``[-2/h**2]``.

The retained computational authority is one canonical, owned ``float64`` CSR
matrix with sorted indices, summed duplicates, and no explicit zeros.  Sparse
and dense inspection return defensive copies.  Application validates finite
real or complex active vectors and returns owned C-contiguous ``float64`` or
``complex128`` arrays.  Wrong semantic types raise :class:`TypeError`; shape,
rank, and finiteness violations raise :class:`ValueError`.

The centered rule targets a second derivative for sufficiently smooth data,
but implementation and execution do not by themselves establish convergence,
physical validation, pedagogical acceptance, or uncertainty quantification.
No physical operator, units, composition, or non-Dirichlet boundary behavior
is introduced.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike
from scipy.sparse import csr_matrix, diags

from projectkoios.physkit.discretization.state_space_1d import (
    HomogeneousDirichletStateSpace1D,
    _validated_numeric_vector,
)

from .discrete_1d import DiscreteLinearOperator1D


class FiniteDifferenceLaplacian1D(DiscreteLinearOperator1D):
    """Immutable centered homogeneous-Dirichlet discrete Laplacian.

    Parameters
    ----------
    state_space : HomogeneousDirichletStateSpace1D
        Immutable domain and codomain.  Its contained grid supplies the
        spacing and its active dimension supplies the matrix shape.

    Attributes
    ----------
    domain : HomogeneousDirichletStateSpace1D
        Exact state-space instance passed to the constructor.
    codomain : HomogeneousDirichletStateSpace1D
        The same exact state-space instance as :attr:`domain`.
    shape : tuple[int, int]
        ``(state_space.dimension, state_space.dimension)``.
    dtype : numpy.dtype
        Exactly ``numpy.dtype(numpy.float64)``.
    matrix : scipy.sparse.csr_matrix
        Defensive copy of the canonical real CSR stencil.

    Raises
    ------
    TypeError
        If ``state_space`` is not a
        :class:`HomogeneousDirichletStateSpace1D`.

    Notes
    -----
    The sign convention is positive second derivative, so the matrix has
    negative diagonal and positive adjacent off-diagonals.  Homogeneous zero
    endpoints are omitted from the active vector.  Domain/codomain
    compatibility is determined by state-space semantic identity rather than
    shape or approximate grid equality.

    Examples
    --------
    >>> from projectkoios.physkit.discretization import (
    ...     HomogeneousDirichletStateSpace1D, UniformGrid1D,
    ... )
    >>> space = HomogeneousDirichletStateSpace1D(UniformGrid1D(0.0, 1.0, 3))
    >>> operator = FiniteDifferenceLaplacian1D(space)
    >>> operator.to_dense()
    array([[-8.]])
    >>> operator @ [1.0]
    array([-8.])

    See Also
    --------
    projectkoios.physkit.discretization.HomogeneousDirichletStateSpace1D
        Boundary and active-vector representation semantics.
    DiscreteLinearOperator1D
        Abstract CSR-matrix operator interface.
    """

    __slots__ = (
        "_state_space",
        "_shape",
        "_dtype",
        "_matrix",
        "_frozen",
    )

    def __init__(self, state_space: HomogeneousDirichletStateSpace1D) -> None:
        if not isinstance(state_space, HomogeneousDirichletStateSpace1D):
            raise TypeError("state_space must be a HomogeneousDirichletStateSpace1D.")

        dimension = state_space.dimension
        main = -2.0 * np.ones(dimension, dtype=np.float64)
        diagonals: list[np.ndarray] = [main]
        offsets = [0]
        # The N=3 case has M=1 and intentionally has no off-diagonal arrays.
        if dimension > 1:
            off = np.ones(dimension - 1, dtype=np.float64)
            diagonals.extend((off, off))
            offsets.extend((-1, 1))
        matrix = diags(
            diagonals,
            offsets,
            shape=(dimension, dimension),
            format="csr",
            dtype=np.float64,
        )
        matrix /= state_space.grid.spacing**2
        matrix.sum_duplicates()
        matrix.sort_indices()
        matrix.eliminate_zeros()

        object.__setattr__(self, "_state_space", state_space)
        object.__setattr__(self, "_shape", (dimension, dimension))
        object.__setattr__(self, "_dtype", np.dtype(np.float64))
        object.__setattr__(self, "_matrix", matrix.copy())
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
    def domain(self) -> HomogeneousDirichletStateSpace1D:
        """Return the exact homogeneous-Dirichlet input state space.

        Returns
        -------
        HomogeneousDirichletStateSpace1D
            Exact constructor argument.
        """
        return self._state_space

    @property
    def codomain(self) -> HomogeneousDirichletStateSpace1D:
        """Return the exact homogeneous-Dirichlet output state space.

        Returns
        -------
        HomogeneousDirichletStateSpace1D
            Same instance as :attr:`domain`.
        """
        return self._state_space

    @property
    def shape(self) -> tuple[int, int]:
        """Return active output-by-input dimensions.

        Returns
        -------
        tuple[int, int]
            Square shape ``(dimension, dimension)``.
        """
        return self._shape

    @property
    def dtype(self) -> np.dtype:
        """Return the real stencil coefficient dtype.

        Returns
        -------
        numpy.dtype
            Exactly ``numpy.dtype(numpy.float64)``.
        """
        return self._dtype

    @property
    def matrix(self) -> csr_matrix:
        """Return a defensive copy of the canonical CSR stencil.

        Returns
        -------
        scipy.sparse.csr_matrix
            New ``float64`` CSR matrix with sorted indices, no duplicate
            entries, and no explicit zeros.  Its buffers do not alias the
            retained computational matrix.
        """
        return self._matrix.copy()

    def apply(self, state: ArrayLike) -> np.ndarray:
        """Apply the centered second-derivative stencil.

        Parameters
        ----------
        state : numpy.typing.ArrayLike
            Finite numeric one-dimensional active vector of shape
            ``(domain.dimension,)``.

        Returns
        -------
        numpy.ndarray
            New owned C-contiguous result.  Real input produces ``float64``;
            complex input produces ``complex128``.  The real stencil acts
            independently on real and imaginary parts.

        Raises
        ------
        TypeError
            If ``state`` is Boolean, object, string, ragged, or otherwise
            nonnumeric.
        ValueError
            If ``state`` has the wrong rank or shape, contains nonfinite
            values, or cannot be finitely represented in its canonical dtype.

        Notes
        -----
        Omitted endpoint terms are the exact homogeneous zero boundary
        values; the input is never reshaped, clipped, or boundary-corrected.
        """
        active = _validated_numeric_vector(state, self._state_space.dimension, "state")
        result = self._matrix @ active
        result_dtype = np.complex128 if np.iscomplexobj(active) else np.float64
        return np.array(result, dtype=result_dtype, order="C", copy=True)
