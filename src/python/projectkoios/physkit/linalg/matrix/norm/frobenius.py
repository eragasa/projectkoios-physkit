"""Evaluate Frobenius norms for canonical PhysKit matrix quantities.

Dense and sparse quantity representations delegate to their exact Optimization owners.
The adapter restores the matrix-entry unit and contains no numerical reduction kernel.
"""

from __future__ import annotations

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.numerics.linalg.numpy.frobenius import NumpyFrobeniusNorm
from projectkoios.numerics.linalg.scipy.frobenius import ScipyCsrFrobeniusNorm

from .base import BaseMatrixQuantityNorm


class FrobeniusMatrixNorm(BaseMatrixQuantityNorm):
    r"""Evaluate the Frobenius norm of one matrix quantity.

    The mathematical identity is

    .. math::

        \lVert A \rVert_F = \sqrt{\sum_{i,j} \lvert A_{ij} \rvert^2}.

    Dense quantities use :class:`NumpyFrobeniusNorm`. Sparse quantities use
    :class:`ScipyCsrFrobeniusNorm` without dense materialization. The inherited
    :meth:`BaseMatrixQuantityNorm.result` method restores the exact entry unit.
    """

    __slots__ = ()

    @property
    def matrix_norm(self) -> MatrixNormEnum:
        """Return the canonical Frobenius semantic identity.

        Returns
        -------
        MatrixNormEnum
            Exact ``MatrixNormEnum.frobenius`` member from Optimization.
        """
        return MatrixNormEnum.frobenius

    def _evaluate_magnitude(self) -> float:
        """Evaluate the unit-free Frobenius magnitude.

        Returns
        -------
        float
            Built-in finite nonnegative Frobenius norm.

        Notes
        -----
        Sparse dispatch obtains one fresh canonical CSR value. It never invokes a dense
        conversion or alternate evaluator.
        """
        if self._matrix_is_sparse():
            return ScipyCsrFrobeniusNorm(self._csr_matrix()).result()
        return NumpyFrobeniusNorm(self._dense_matrix()).result()


__all__ = ["FrobeniusMatrixNorm"]
