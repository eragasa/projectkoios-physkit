"""Evaluate maximum absolute entries for canonical matrix quantities.

Dense and sparse quantity representations delegate to their exact Optimization owners.
The adapter restores the matrix-entry unit and contains no local reduction kernel.
"""

from __future__ import annotations

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.numerics.linalg.numpy.maximum_absolute_entry import (
    NumpyMaximumAbsoluteEntryNorm,
)
from projectkoios.numerics.linalg.scipy.maximum_absolute_entry import (
    ScipyCsrMaximumAbsoluteEntryNorm,
)

from .base import BaseMatrixQuantityNorm


class MaximumAbsoluteEntryMatrixNorm(BaseMatrixQuantityNorm):
    r"""Evaluate the maximum absolute entry of one matrix quantity.

    The mathematical identity is

    .. math::

        \lVert A \rVert_{\max} = \max_{i,j} \lvert A_{ij} \rvert.

    Dense quantities use :class:`NumpyMaximumAbsoluteEntryNorm`. Sparse quantities use
    :class:`ScipyCsrMaximumAbsoluteEntryNorm`, which includes represented implicit zeros
    without dense materialization. The inherited result boundary restores the exact
    entry unit.
    """

    __slots__ = ()

    @property
    def matrix_norm(self) -> MatrixNormEnum:
        """Return the canonical maximum-entry semantic identity.

        Returns
        -------
        MatrixNormEnum
            Exact ``MatrixNormEnum.maximum_absolute_entry`` member from Optimization.
        """
        return MatrixNormEnum.maximum_absolute_entry

    def _evaluate_magnitude(self) -> float:
        """Evaluate the unit-free maximum absolute entry.

        Returns
        -------
        float
            Built-in finite nonnegative maximum entry magnitude.

        Notes
        -----
        Sparse dispatch obtains one fresh canonical CSR value. It never invokes a dense
        conversion or alternate evaluator.
        """
        if self._matrix_is_sparse():
            return ScipyCsrMaximumAbsoluteEntryNorm(self._csr_matrix()).result()
        return NumpyMaximumAbsoluteEntryNorm(self._dense_matrix()).result()


__all__ = ["MaximumAbsoluteEntryMatrixNorm"]
