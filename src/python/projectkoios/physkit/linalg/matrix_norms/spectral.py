"""Spectral norms of canonical PhysKit matrix quantities."""

from __future__ import annotations

from enum import Enum, auto

import numpy as np

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.numerics.linalg.numpy.spectral_norm import NumpySpectralNorm
from projectkoios.numerics.linalg.scipy.spectral_norm import (
    ScipyCsrSpectralNorm,
    ScipyCsrSpectralNormControls,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    SparseMatrixQuantity,
)

from .base import BaseMatrixQuantityNorm


class SparseSpectralNormOptionsEnum(Enum):
    """Select one explicit spectral path for a sparse matrix quantity."""

    complete = auto()
    iterative = auto()


class SpectralMatrixNorm(BaseMatrixQuantityNorm):
    """Evaluate one matrix quantity's largest singular value."""

    __slots__ = ("_sparse_spectral_controls", "_sparse_spectral_option")

    def __init__(
        self,
        matrix: (
            MatrixQuantity
            | ComplexMatrixQuantity
            | SparseMatrixQuantity
            | ComplexSparseMatrixQuantity
        ),
        /,
        *,
        sparse_spectral_option: SparseSpectralNormOptionsEnum | None = None,
        sparse_spectral_controls: ScipyCsrSpectralNormControls | None = None,
    ) -> None:
        """Retain one matrix quantity and validate sparse spectral selection."""
        super().__init__(matrix)
        if (
            sparse_spectral_option is not None
            and type(sparse_spectral_option) is not SparseSpectralNormOptionsEnum
        ):
            raise TypeError(
                "sparse_spectral_option must be SparseSpectralNormOptionsEnum or None"
            )
        if (
            sparse_spectral_controls is not None
            and type(sparse_spectral_controls) is not ScipyCsrSpectralNormControls
        ):
            raise TypeError(
                "sparse_spectral_controls must be ScipyCsrSpectralNormControls or None"
            )

        if not self._matrix_is_sparse():
            if sparse_spectral_option is not None:
                raise ValueError(
                    "dense spectral evaluation does not accept a sparse option"
                )
            if sparse_spectral_controls is not None:
                raise ValueError(
                    "dense spectral evaluation does not accept sparse controls"
                )
        elif sparse_spectral_option is None:
            raise ValueError("sparse spectral evaluation requires an explicit option")
        elif sparse_spectral_option is SparseSpectralNormOptionsEnum.complete:
            if sparse_spectral_controls is not None:
                raise ValueError(
                    "complete sparse spectral evaluation does not accept iterative "
                    "controls"
                )
        elif sparse_spectral_controls is None:
            raise ValueError(
                "iterative sparse spectral evaluation requires explicit controls"
            )

        self._sparse_spectral_option = sparse_spectral_option
        self._sparse_spectral_controls = sparse_spectral_controls

    @property
    def matrix_norm(self) -> MatrixNormEnum:
        """Return the canonical spectral semantic identity."""
        return MatrixNormEnum.spectral

    @property
    def sparse_spectral_option(self) -> SparseSpectralNormOptionsEnum | None:
        """Return the retained sparse spectral selection, when applicable."""
        return self._sparse_spectral_option

    @property
    def sparse_spectral_controls(self) -> ScipyCsrSpectralNormControls | None:
        """Return the retained iterative controls, when applicable."""
        return self._sparse_spectral_controls

    def _evaluate_magnitude(self) -> float:
        """Delegate to the selected complete or iterative numerical owner."""
        if not self._matrix_is_sparse():
            return NumpySpectralNorm(self._dense_matrix()).result()
        if self._sparse_spectral_option is SparseSpectralNormOptionsEnum.complete:
            return self._complete_sparse_result()
        return self._iterative_sparse_result()

    def _complete_sparse_result(self) -> float:
        """Materialize CSR explicitly and evaluate a complete dense SVD."""
        matrix = self._csr_matrix()
        row_count, column_count = matrix.shape
        if row_count == 0 or column_count == 0:
            raise ValueError("matrix must be nonempty")
        with np.errstate(over="ignore", invalid="ignore"):
            dense = matrix.toarray()
        if not bool(np.isfinite(dense).all()):
            raise FloatingPointError("CSR materialization produced nonfinite entries")
        return NumpySpectralNorm(dense).result()

    def _iterative_sparse_result(self) -> float:
        """Evaluate CSR without dense conversion or fallback."""
        controls = self._sparse_spectral_controls
        if controls is None:
            raise RuntimeError("iterative sparse spectral controls are unavailable")
        return ScipyCsrSpectralNorm(self._csr_matrix(), controls).result()


__all__ = [
    "ScipyCsrSpectralNormControls",
    "SparseSpectralNormOptionsEnum",
    "SpectralMatrixNorm",
]
