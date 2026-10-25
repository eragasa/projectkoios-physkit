"""Define quantity-aware matrix-norm contracts and semantic dispatch.

The module owns PhysKit's nominal quantity boundary and exact unit restoration.
Optimization supplies the canonical norm enum and every unit-free numerical evaluator.
"""

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
    # Runtime import would create a cycle because spectral evaluators derive from the
    # base declared in this module. Postponed annotations retain the public signature.
    from .spectral import SparseSpectralNormOptionsEnum


class BaseMatrixQuantityNorm(ABC):
    """Define shared behavior for one norm over a canonical matrix quantity.

    Parameters
    ----------
    matrix
        Exact real or complex dense or sparse PhysKit matrix quantity. The evaluator
        retains this object by identity.

    Notes
    -----
    The class owns quantity representation access and unit restoration. Concrete
    subclasses invoke Optimization evaluators for the unit-free magnitude. Raw NumPy
    arrays and SciPy sparse objects are not accepted at this boundary.
    """

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
        """Retain one exact dense or sparse matrix quantity.

        Parameters
        ----------
        matrix
            Canonical real or complex dense or sparse matrix quantity.

        Raises
        ------
        TypeError
            If ``matrix`` is not one of the four nominal matrix quantity classes.
        """
        # Nominal ownership is intentional. Structurally similar values do not carry
        # PhysKit's immutable storage and exact unit contracts.
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
        """Return the exact retained matrix quantity.

        Returns
        -------
        MatrixQuantity | ComplexMatrixQuantity | SparseMatrixQuantity |
        ComplexSparseMatrixQuantity
            Matrix quantity supplied at construction, preserved by identity.
        """
        return self._matrix

    @property
    @abstractmethod
    def matrix_norm(self) -> MatrixNormEnum:
        """Return the exact Optimization norm semantic identity.

        Returns
        -------
        MatrixNormEnum
            Canonical enum member implemented by the concrete evaluator.
        """
        raise NotImplementedError

    def result(self) -> ScalarQuantity:
        """Evaluate the norm and restore the exact matrix-entry unit.

        Returns
        -------
        ScalarQuantity
            Built-in finite nonnegative owner result paired with the exact unit object
            retained by the matrix quantity.

        Notes
        -----
        Evaluation executes on every call. The base performs no caching, unit
        conversion, numerical reduction, retry, or fallback.
        """
        return ScalarQuantity(self._evaluate_magnitude(), self._matrix.unit)

    def _matrix_is_sparse(self) -> bool:
        """Return whether the retained quantity owns a CSR representation.

        Returns
        -------
        bool
            ``True`` for real or complex sparse matrix quantities and ``False`` for
            real or complex dense matrix quantities.
        """
        return isinstance(
            self._matrix,
            SparseMatrixQuantity | ComplexSparseMatrixQuantity,
        )

    def _dense_matrix(self) -> DenseMatrix:
        """Return the exact immutable magnitude of a dense quantity.

        Returns
        -------
        DenseMatrix
            Exact binary64 or complex128 ndarray owned by the retained quantity.

        Raises
        ------
        TypeError
            If the retained quantity uses sparse storage.
        """
        if isinstance(self._matrix, MatrixQuantity | ComplexMatrixQuantity):
            return self._matrix.magnitude
        raise TypeError("matrix must be a dense matrix quantity")

    def _csr_matrix(self) -> CsrMatrix:
        """Return a fresh canonical CSR value from a sparse quantity.

        Returns
        -------
        CsrMatrix
            Fresh binary64 or complex128 CSR array defined by the quantity's explicit
            copy boundary.

        Raises
        ------
        TypeError
            If the retained quantity uses dense storage.
        """
        if isinstance(
            self._matrix,
            SparseMatrixQuantity | ComplexSparseMatrixQuantity,
        ):
            return self._matrix.to_csr()
        raise TypeError("matrix must be a sparse matrix quantity")

    @abstractmethod
    def _evaluate_magnitude(self) -> float:
        """Return one unit-free magnitude from the exact numerical owner.

        Returns
        -------
        float
            Built-in finite nonnegative matrix-norm magnitude.
        """
        raise NotImplementedError


@dataclass(frozen=True, slots=True, eq=False)
class MatrixNorm:
    """Dispatch one matrix quantity to a selected quantity-aware norm.

    Parameters
    ----------
    matrix
        Exact real or complex dense or sparse PhysKit matrix quantity.
    matrix_norm
        Exact canonical Optimization norm enum member.
    sparse_spectral_option
        Exact sparse spectral path selector. Supply only for spectral evaluation of a
        sparse matrix quantity.
    sparse_spectral_controls
        Exact Optimization controls. Supply only with the ``iterative`` sparse spectral
        option.

    Notes
    -----
    Sparse spectral arguments are keyword-only. The facade is frozen, slotted, and
    identity-based. It delegates one call without numerical arithmetic, unit conversion,
    dense conversion, result caching, retry, fallback, or exception translation.
    """

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
        """Validate public identities and install immutable dispatch state.

        Raises
        ------
        TypeError
            If the matrix, enum, sparse option, or controls object has an unsupported
            nominal identity.
        ValueError
            If a nonspectral norm receives a sparse spectral option or controls.
        """
        # Local imports break the base-to-concrete import cycle while keeping concrete
        # implementation identities out of module-level mutable registries.
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

        # Exact enum and control types prevent lookalike values from silently selecting
        # a numerical algorithm owned by another package.
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

        # Mapping proxies preserve an explicit closed dispatch table. Underscore names
        # and repr/compare exclusion keep these implementation details non-public.
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

        # Sparse policy has meaning only for the spectral implementation. Rejecting it
        # here prevents silent option loss in Frobenius and maximum-entry evaluation.
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
        """Evaluate the selected norm and return its unit-aware scalar.

        Returns
        -------
        ScalarQuantity
            Norm magnitude paired with the exact unit object retained by ``matrix``.

        Raises
        ------
        TypeError
            If concrete-owner representation validation fails.
        ValueError
            If spectral sparse arguments or matrix shape violate the selected path.
        FloatingPointError
            If materialization or owner arithmetic produces a nonfinite result.
        RuntimeError
            If the selected dense or iterative decomposition does not converge.
        """
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
