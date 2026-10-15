"""Abstract finite-matrix operators on one-dimensional state spaces.

:class:`DiscreteLinearOperator1D` specializes
:class:`~projectkoios.physkit.operators.base.LinearOperator` by requiring a
canonical SciPy CSR matrix. Dense inspection is derived from the public defensive
CSR copy, not from an independent representation, and returns owned C-contiguous data.
Concrete subclasses define semantic domain/codomain validation, matrix sign,
shape, and dtype.

This module introduces no grid, boundary condition, stencil, composition, or
physical model.  Matrix inspection alone provides neither numerical nor
physical validation, pedagogical acceptance, nor uncertainty quantification.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from scipy.sparse import csr_matrix

from .base import LinearOperator


class DiscreteLinearOperator1D(LinearOperator, ABC):
    """Abstract one-dimensional operator with a canonical CSR matrix.

    Attributes
    ----------
    matrix : scipy.sparse.csr_matrix
        Defensive copy of the concrete operator's canonical matrix.

    Notes
    -----
    Domain and codomain semantics, application validation, sign, shape, and
    dtype remain responsibilities of concrete subclasses.  The matrix is the
    computational representation; :meth:`to_dense` is inspection only.

    See Also
    --------
    LinearOperator
        General semantic operator interface.
    projectkoios.physkit.operators.FiniteDifferenceLaplacian1D
        Centered homogeneous-Dirichlet concrete implementation.
    """

    @property
    @abstractmethod
    def matrix(self) -> csr_matrix:
        """Return a defensive copy of the canonical sparse matrix.

        Returns
        -------
        scipy.sparse.csr_matrix
            New CSR matrix whose mutable buffers do not alias retained
            operator storage.
        """
        raise NotImplementedError

    def to_dense(self) -> np.ndarray:
        """Return an owned dense matrix derived from CSR.

        Returns
        -------
        numpy.ndarray
            New C-contiguous dense array with this operator's matrix dtype
            and shape.

        Notes
        -----
        The dense form is not retained and is never a second computational
        authority.
        """
        return np.array(self.matrix.toarray(), order="C", copy=True)
