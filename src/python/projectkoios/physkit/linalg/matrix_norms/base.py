"""Common contracts and semantic dispatch for matrix-quantity norms."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.numerics.linalg.scipy.spectral_norm import (
    ScipyCsrSpectralNormControls,
)
from projectkoios.numerics.numpy.typing import DenseMatrix
from projectkoios.numerics.scipy.typing import CsrMatrix
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    ScalarQuantity,
    SparseMatrixQuantity,
)

if TYPE_CHECKING:
    from .spectral import SparseSpectralNormOptionsEnum


class BaseMatrixQuantityNorm(ABC):
    """Define one norm over a canonical PhysKit matrix quantity."""

    __slots__ = ("_matrix",)

    def __init__(
        self,
        matrix: (
            MatrixQuantity
            | ComplexMatrixQuantity
            | SparseMatrixQuantity
            | ComplexSparseMatrixQuantity
        ),
        /,
    ) -> None:
        """Retain one exact dense or sparse matrix quantity."""
        if not isinstance(
            matrix,
            MatrixQuantity
            | ComplexMatrixQuantity
            | SparseMatrixQuantity
            | ComplexSparseMatrixQuantity,
        ):
            raise TypeError("matrix must be a real or complex matrix quantity")
        self._matrix = matrix

    @property
    def matrix(
        self,
    ) -> (
        MatrixQuantity
        | ComplexMatrixQuantity
        | SparseMatrixQuantity
        | ComplexSparseMatrixQuantity
    ):
        """Return the exact retained matrix quantity."""
        return self._matrix

    @property
    @abstractmethod
    def matrix_norm(self) -> MatrixNormEnum:
        """Return the exact Optimization norm semantic identity."""
        raise NotImplementedError

    def result(self) -> ScalarQuantity:
        """Return the norm magnitude with the exact retained matrix unit."""
        return ScalarQuantity(self._evaluate_magnitude(), self._matrix.unit)

    def _matrix_is_sparse(self) -> bool:
        """Return whether the retained quantity owns a CSR representation."""
        return isinstance(
            self._matrix,
            SparseMatrixQuantity | ComplexSparseMatrixQuantity,
        )

    def _dense_matrix(self) -> DenseMatrix:
        """Return the exact dense magnitude from a dense matrix quantity."""
        if isinstance(self._matrix, MatrixQuantity | ComplexMatrixQuantity):
            return self._matrix.magnitude
        raise TypeError("matrix must be a dense matrix quantity")

    def _csr_matrix(self) -> CsrMatrix:
        """Return a fresh canonical CSR value from a sparse matrix quantity."""
        if isinstance(
            self._matrix,
            SparseMatrixQuantity | ComplexSparseMatrixQuantity,
        ):
            return self._matrix.to_csr()
        raise TypeError("matrix must be a sparse matrix quantity")

    @abstractmethod
    def _evaluate_magnitude(self) -> float:
        """Return one unit-free magnitude from the exact numerical owner."""
        raise NotImplementedError


@dataclass(frozen=True, slots=True, eq=False)
class MatrixNorm:
    """Dispatch one matrix quantity to its selected quantity-aware norm."""

    matrix: (
        MatrixQuantity
        | ComplexMatrixQuantity
        | SparseMatrixQuantity
        | ComplexSparseMatrixQuantity
    )
    matrix_norm: MatrixNormEnum
    _: KW_ONLY
    sparse_spectral_option: SparseSpectralNormOptionsEnum | None = None
    sparse_spectral_controls: ScipyCsrSpectralNormControls | None = None
    _norm_names: Mapping[MatrixNormEnum, str] = field(
        init=False,
        repr=False,
        compare=False,
    )
    _implementations: Mapping[MatrixNormEnum, type[BaseMatrixQuantityNorm]] = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Install immutable dispatch state and validate facade-owned choices."""
        from .frobenius import FrobeniusMatrixNorm
        from .maximum_absolute_entry import MaximumAbsoluteEntryMatrixNorm
        from .spectral import SparseSpectralNormOptionsEnum, SpectralMatrixNorm

        if not isinstance(
            self.matrix,
            MatrixQuantity
            | ComplexMatrixQuantity
            | SparseMatrixQuantity
            | ComplexSparseMatrixQuantity,
        ):
            raise TypeError("matrix must be a real or complex matrix quantity")
        if type(self.matrix_norm) is not MatrixNormEnum:
            raise TypeError("matrix_norm must be MatrixNormEnum")
        if (
            self.sparse_spectral_option is not None
            and type(self.sparse_spectral_option) is not SparseSpectralNormOptionsEnum
        ):
            raise TypeError(
                "sparse_spectral_option must be SparseSpectralNormOptionsEnum or None"
            )
        if (
            self.sparse_spectral_controls is not None
            and type(self.sparse_spectral_controls) is not ScipyCsrSpectralNormControls
        ):
            raise TypeError(
                "sparse_spectral_controls must be ScipyCsrSpectralNormControls or None"
            )

        object.__setattr__(
            self,
            "_norm_names",
            MappingProxyType(
                {
                    MatrixNormEnum.frobenius: "frobenius",
                    MatrixNormEnum.spectral: "spectral",
                    MatrixNormEnum.maximum_absolute_entry: "maximum-absolute-entry",
                }
            ),
        )
        object.__setattr__(
            self,
            "_implementations",
            MappingProxyType(
                {
                    MatrixNormEnum.frobenius: FrobeniusMatrixNorm,
                    MatrixNormEnum.spectral: SpectralMatrixNorm,
                    MatrixNormEnum.maximum_absolute_entry: (
                        MaximumAbsoluteEntryMatrixNorm
                    ),
                }
            ),
        )

        if self.matrix_norm is not MatrixNormEnum.spectral:
            norm_name = self._norm_names[self.matrix_norm]
            if self.sparse_spectral_option is not None:
                raise ValueError(
                    f"{norm_name} evaluation does not accept a sparse spectral option"
                )
            if self.sparse_spectral_controls is not None:
                raise ValueError(
                    f"{norm_name} evaluation does not accept sparse spectral controls"
                )

    def result(self) -> ScalarQuantity:
        """Return the selected norm with the exact retained matrix unit."""
        implementation = self._implementations[self.matrix_norm]
        if self.matrix_norm is MatrixNormEnum.spectral:
            evaluation = implementation(
                self.matrix,
                sparse_spectral_option=self.sparse_spectral_option,
                sparse_spectral_controls=self.sparse_spectral_controls,
            )
        else:
            evaluation = implementation(self.matrix)
        return evaluation.result()


__all__ = ["BaseMatrixQuantityNorm", "MatrixNorm", "MatrixNormEnum"]
