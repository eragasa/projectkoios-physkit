"""Maximum-absolute-entry norms of canonical PhysKit matrix quantities."""

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
    """Evaluate one matrix quantity's maximum absolute entry."""

    __slots__ = ()

    @property
    def matrix_norm(self) -> MatrixNormEnum:
        """Return the canonical maximum-absolute-entry semantic identity."""
        return MatrixNormEnum.maximum_absolute_entry

    def _evaluate_magnitude(self) -> float:
        """Delegate the represented magnitude to its exact numerical owner."""
        if self._matrix_is_sparse():
            return ScipyCsrMaximumAbsoluteEntryNorm(self._csr_matrix()).result()
        return NumpyMaximumAbsoluteEntryNorm(self._dense_matrix()).result()


__all__ = ["MaximumAbsoluteEntryMatrixNorm"]
