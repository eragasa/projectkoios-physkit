"""Evaluate spectral norms for canonical PhysKit matrix quantities.

The module owns explicit dense-versus-sparse path selection and exact unit restoration.
Optimization owns complete singular-value decomposition, sparse direct arithmetic,
ARPACK iteration, residual qualification, and numerical failure classification.
"""

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
    """Select one explicit spectral path for a sparse matrix quantity.

    Attributes
    ----------
    complete
        Materialize one canonical CSR value exactly once and invoke the complete dense
        singular-value decomposition owner. Iterative controls are invalid.
    iterative
        Invoke the sparse Optimization owner with exact explicit controls. Dense
        materialization and complete-SVD fallback are prohibited.

    Notes
    -----
    The enum values are opaque and are not strings. Sparse spectral evaluation has no
    implicit or default path.
    """

    complete = auto()
    iterative = auto()


class SpectralMatrixNorm(BaseMatrixQuantityNorm):
    r"""Evaluate the spectral norm of one matrix quantity.

    The mathematical identity is

    .. math::

        \lVert A \rVert_2 = \sigma_{\max}(A).

    Parameters
    ----------
    matrix
        Exact real or complex dense or sparse PhysKit matrix quantity.
    sparse_spectral_option
        ``None`` for a dense quantity. A sparse quantity requires exact ``complete`` or
        ``iterative`` selection.
    sparse_spectral_controls
        ``None`` for dense or complete evaluation. Iterative sparse evaluation requires
        one exact :class:`ScipyCsrSpectralNormControls` instance.

    Notes
    -----
    Dense quantities use complete NumPy singular-value decomposition. Complete sparse
    evaluation contains the package's sole CSR-to-dense boundary. Iterative sparse
    evaluation remains CSR and never retries through the complete path.
    """

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
        """Retain a matrix quantity and validate spectral path selection.

        Parameters
        ----------
        matrix
            Exact real or complex dense or sparse matrix quantity.
        sparse_spectral_option
            Exact sparse path selector or ``None`` for dense evaluation.
        sparse_spectral_controls
            Exact iterative controls or ``None`` for dense and complete evaluation.

        Raises
        ------
        TypeError
            If the matrix, option, or controls object has an unsupported nominal
            identity.
        ValueError
            If sparse arguments conflict with the matrix representation or selected
            path.
        """
        super().__init__(matrix)

        # Exact identities prevent strings or lookalike control records from selecting
        # a numerical algorithm implicitly.
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

        # The decision table rejects every ignored or inferred sparse policy. Dense
        # evaluation has one complete path; sparse evaluation requires an explicit path.
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
        """Return the canonical spectral semantic identity.

        Returns
        -------
        MatrixNormEnum
            Exact ``MatrixNormEnum.spectral`` member from Optimization.
        """
        return MatrixNormEnum.spectral

    @property
    def sparse_spectral_option(self) -> SparseSpectralNormOptionsEnum | None:
        """Return the exact retained sparse path selector.

        Returns
        -------
        SparseSpectralNormOptionsEnum | None
            ``None`` for dense evaluation; exact ``complete`` or ``iterative`` member
            for sparse evaluation.
        """
        return self._sparse_spectral_option

    @property
    def sparse_spectral_controls(self) -> ScipyCsrSpectralNormControls | None:
        """Return the exact retained iterative controls.

        Returns
        -------
        ScipyCsrSpectralNormControls | None
            Exact controls for iterative sparse evaluation; ``None`` for dense and
            complete sparse evaluation.
        """
        return self._sparse_spectral_controls

    def _evaluate_magnitude(self) -> float:
        """Evaluate the unit-free spectral magnitude through the selected owner.

        Returns
        -------
        float
            Built-in finite nonnegative largest singular value.

        Raises
        ------
        ValueError
            If complete sparse evaluation receives an empty matrix.
        FloatingPointError
            If complete materialization or owner arithmetic becomes nonfinite.
        RuntimeError
            If the complete or iterative numerical owner does not converge.
        """
        if not self._matrix_is_sparse():
            return NumpySpectralNorm(self._dense_matrix()).result()
        if self._sparse_spectral_option is SparseSpectralNormOptionsEnum.complete:
            return self._complete_sparse_result()
        return self._iterative_sparse_result()

    def _complete_sparse_result(self) -> float:
        """Materialize CSR once and invoke complete dense SVD.

        Returns
        -------
        float
            Built-in finite nonnegative largest singular value.

        Raises
        ------
        ValueError
            If either matrix axis has length zero.
        FloatingPointError
            If explicit materialization produces a nonfinite entry.
        RuntimeError
            If the dense singular-value decomposition does not converge.
        """
        matrix = self._csr_matrix()
        row_count, column_count = matrix.shape
        if row_count == 0 or column_count == 0:
            # Reject emptiness before allocating the authorized dense representation.
            raise ValueError("matrix must be nonempty")

        # Complete is the only path that authorizes toarray(). Floating warnings from
        # representability are classified by the explicit finite check below.
        with np.errstate(over="ignore", invalid="ignore"):
            dense = matrix.toarray()
        if not bool(np.isfinite(dense).all()):
            raise FloatingPointError("CSR materialization produced nonfinite entries")
        return NumpySpectralNorm(dense).result()

    def _iterative_sparse_result(self) -> float:
        """Invoke the CSR owner without dense conversion or fallback.

        Returns
        -------
        float
            Built-in finite nonnegative largest singular value.

        Raises
        ------
        RuntimeError
            If validated construction state lacks controls or the sparse numerical owner
            does not converge.
        FloatingPointError
            If sparse owner arithmetic or residual qualification becomes nonfinite.
        """
        controls = self._sparse_spectral_controls
        if controls is None:
            # The public constructor rejects this state. The guard prevents a missing
            # control record from reaching the numerical owner after constructor
            # execution has been bypassed.
            raise RuntimeError("iterative sparse spectral controls are unavailable")
        return ScipyCsrSpectralNorm(self._csr_matrix(), controls).result()


__all__ = [
    "ScipyCsrSpectralNormControls",
    "SparseSpectralNormOptionsEnum",
    "SpectralMatrixNorm",
]
