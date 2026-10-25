"""Frobenius norms of canonical PhysKit matrix quantities."""

from __future__ import annotations

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.numerics.linalg.numpy.frobenius import NumpyFrobeniusNorm
from projectkoios.numerics.linalg.scipy.frobenius import ScipyCsrFrobeniusNorm

from .base import BaseMatrixQuantityNorm


class FrobeniusMatrixNorm(BaseMatrixQuantityNorm):
    """Evaluate one matrix quantity's Frobenius norm."""

    __slots__ = ()

    @property
    def matrix_norm(self) -> MatrixNormEnum:
        """Return the canonical Frobenius semantic identity."""
        return MatrixNormEnum.frobenius

    def _evaluate_magnitude(self) -> float:
        """Delegate the represented magnitude to its exact numerical owner."""
        if self._matrix_is_sparse():
            return ScipyCsrFrobeniusNorm(self._csr_matrix()).result()
        return NumpyFrobeniusNorm(self._dense_matrix()).result()


__all__ = ["FrobeniusMatrixNorm"]
