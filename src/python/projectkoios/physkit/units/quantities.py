"""Immutable typed quantities for public represented operators and analyses.

Pint owns unit parsing, dimensional compatibility, and conversion. Project records own
strict scalar and array types, operational immutability, and the distinction between a
physical unit and the first-class :class:`Unitless` unit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pint
from scipy import sparse  # type: ignore[import-untyped]

from projectkoios.physkit.numerics.typing.numpy.arrays import (
    ComplexMatrix,
    ComplexVector,
    IntegerVector,
    RealMatrix,
    RealVector,
)
from projectkoios.physkit.numerics.typing.scipy.sparse import SparseMatrix

__all__ = [
    "MODEL_SYSTEM_PINT_REGISTRY",
    "ModelSystemUnitConverter",
    "ComplexMatrixQuantity",
    "ComplexSparseMatrixQuantity",
    "ComplexVectorQuantity",
    "MatrixQuantity",
    "ModelSystemQuantity",
    "ModelSystemUnit",
    "PhysicalUnit",
    "PintUnitConverter",
    "ScalarQuantity",
    "SparseMatrixQuantity",
    "Unitless",
    "VectorQuantity",
]

MODEL_SYSTEM_PINT_REGISTRY = pint.UnitRegistry()
"""Fixed Pint registry used by public model-system unit records and actions."""


class ModelSystemUnit:
    """Nominal runtime base for physical and explicitly unitless values."""

    if TYPE_CHECKING:

        @property
        def expression(self) -> str:
            """Return the read-only unit expression shared by concrete units."""
            ...


@dataclass(frozen=True, slots=True)
class PhysicalUnit(ModelSystemUnit):
    """Identify one non-dimensionless Pint unit expression.

    Parameters
    ----------
    expression
        Nonempty Pint unit expression. Dimensional validity is checked by
        :class:`PintUnitConverter` at the conversion boundary.
    """

    expression: str

    def __post_init__(self) -> None:
        if not isinstance(self.expression, str):
            raise TypeError("expression must be a string")
        if not self.expression:
            raise ValueError("expression must be nonempty")
        try:
            parsed = MODEL_SYSTEM_PINT_REGISTRY.Unit(self.expression)
        except (
            pint.errors.DefinitionSyntaxError,
            pint.errors.UndefinedUnitError,
        ) as error:
            raise ValueError("expression must identify a valid Pint unit") from error
        if parsed.dimensionless:
            raise ValueError("dimensionless values require Unitless")


@dataclass(frozen=True, slots=True)
class Unitless(ModelSystemUnit):
    """Represent the first-class unit of a genuinely dimensionless quantity."""

    @property
    def expression(self) -> str:
        """Return Pint's canonical dimensionless unit expression."""
        return "dimensionless"


class ModelSystemQuantity:
    """Nominal runtime base for immutable unit-bearing numerical values."""


@dataclass(frozen=True, slots=True)
class ScalarQuantity(ModelSystemQuantity):
    """Represent one finite binary64 scalar and its explicit unit."""

    magnitude: float
    unit: ModelSystemUnit

    def __post_init__(self) -> None:
        if type(self.magnitude) is not float:
            raise TypeError("magnitude must be a built-in float")
        if not np.isfinite(self.magnitude):
            raise ValueError("magnitude must be finite")
        if not isinstance(self.unit, PhysicalUnit | Unitless):
            raise TypeError("unit must be PhysicalUnit or Unitless")


@dataclass(frozen=True, slots=True, eq=False)
class VectorQuantity(ModelSystemQuantity):
    """Represent one immutable finite binary64 vector and its explicit unit."""

    magnitude: RealVector
    unit: ModelSystemUnit

    def __post_init__(self) -> None:
        if not isinstance(self.magnitude, np.ndarray):
            raise TypeError("magnitude must be a numpy.ndarray")
        if self.magnitude.dtype.kind not in "iuf":
            raise TypeError(
                "magnitude must contain real numeric values excluding booleans"
            )
        values = np.asarray(self.magnitude, dtype=np.float64)
        if values.ndim != 1 or not np.all(np.isfinite(values)):
            raise ValueError("magnitude must be a finite vector")
        if not isinstance(self.unit, PhysicalUnit | Unitless):
            raise TypeError("unit must be PhysicalUnit or Unitless")
        immutable = np.frombuffer(values.astype("<f8").tobytes(), dtype="<f8")
        object.__setattr__(self, "magnitude", immutable)


@dataclass(frozen=True, slots=True, eq=False)
class ComplexVectorQuantity(ModelSystemQuantity):
    """Represent one immutable finite complex128 vector and its explicit unit."""

    magnitude: ComplexVector
    unit: ModelSystemUnit

    def __post_init__(self) -> None:
        if not isinstance(self.magnitude, np.ndarray):
            raise TypeError("magnitude must be a numpy.ndarray")
        if self.magnitude.dtype.kind not in "iufc":
            raise TypeError("magnitude must contain numeric values excluding booleans")
        values = np.asarray(self.magnitude, dtype=np.complex128)
        if values.ndim != 1:
            raise ValueError("magnitude must be a vector")
        if not np.all(np.isfinite(values.real)) or not np.all(np.isfinite(values.imag)):
            raise ValueError("magnitude must contain finite complex values")
        if not isinstance(self.unit, PhysicalUnit | Unitless):
            raise TypeError("unit must be PhysicalUnit or Unitless")
        immutable = np.frombuffer(values.astype("<c16").tobytes(), dtype="<c16")
        object.__setattr__(self, "magnitude", immutable)


@dataclass(frozen=True, slots=True, eq=False)
class MatrixQuantity(ModelSystemQuantity):
    """Represent one immutable finite binary64 matrix and its explicit unit."""

    magnitude: RealMatrix
    unit: ModelSystemUnit

    def __post_init__(self) -> None:
        if not isinstance(self.magnitude, np.ndarray):
            raise TypeError("magnitude must be a numpy.ndarray")
        if self.magnitude.dtype.kind not in "iuf":
            raise TypeError(
                "magnitude must contain real numeric values excluding booleans"
            )
        values = np.asarray(self.magnitude, dtype=np.float64)
        if values.ndim != 2 or not np.all(np.isfinite(values)):
            raise ValueError("magnitude must be a finite matrix")
        if not isinstance(self.unit, PhysicalUnit | Unitless):
            raise TypeError("unit must be PhysicalUnit or Unitless")
        immutable = np.frombuffer(values.astype("<f8").tobytes(), dtype="<f8")
        object.__setattr__(self, "magnitude", immutable.reshape(values.shape))


@dataclass(frozen=True, slots=True, eq=False)
class ComplexMatrixQuantity(ModelSystemQuantity):
    """Represent one immutable finite complex128 matrix and its explicit unit.

    Numeric NumPy arrays with integer, floating, or complex dtype are accepted.
    Boolean and nonnumeric arrays are rejected rather than coerced. Stored values are
    copied into immutable byte-backed complex128 storage.
    """

    magnitude: ComplexMatrix
    unit: ModelSystemUnit

    def __post_init__(self) -> None:
        """Validate finite numeric matrix input and install immutable storage."""
        if not isinstance(self.magnitude, np.ndarray):
            raise TypeError("magnitude must be a numpy.ndarray")
        if self.magnitude.dtype.kind not in "iufc":
            raise TypeError("magnitude must contain numeric values excluding booleans")
        values = np.asarray(self.magnitude, dtype=np.complex128)
        if values.ndim != 2:
            raise ValueError("magnitude must be a matrix")
        if not np.all(np.isfinite(values.real)) or not np.all(np.isfinite(values.imag)):
            raise ValueError("magnitude must contain finite complex values")
        if not isinstance(self.unit, PhysicalUnit | Unitless):
            raise TypeError("unit must be PhysicalUnit or Unitless")
        immutable = np.frombuffer(values.astype("<c16").tobytes(), dtype="<c16")
        object.__setattr__(self, "magnitude", immutable.reshape(values.shape))


@dataclass(frozen=True, slots=True, eq=False)
class SparseMatrixQuantity(ModelSystemQuantity):
    """Represent one immutable finite binary64 matrix in canonical CSR storage.

    The three CSR component arrays are copied into immutable byte-backed NumPy
    arrays. :meth:`to_csr` returns a fresh SciPy sparse array, so callers cannot
    mutate the maintained quantity through the returned representation.
    """

    data: RealVector
    column_indices: IntegerVector
    row_offsets: IntegerVector
    shape: tuple[int, int]
    unit: ModelSystemUnit

    def __post_init__(self) -> None:
        if not isinstance(self.data, np.ndarray):
            raise TypeError("data must be a numpy.ndarray")
        if not isinstance(self.column_indices, np.ndarray):
            raise TypeError("column_indices must be a numpy.ndarray")
        if not isinstance(self.row_offsets, np.ndarray):
            raise TypeError("row_offsets must be a numpy.ndarray")
        if self.data.dtype.kind not in "iuf":
            raise TypeError("data must contain real numeric values excluding booleans")
        if self.column_indices.dtype.kind not in "iu":
            raise TypeError("column_indices must contain integers excluding booleans")
        if self.row_offsets.dtype.kind not in "iu":
            raise TypeError("row_offsets must contain integers excluding booleans")
        data = np.asarray(self.data, dtype=np.float64)
        columns = np.asarray(self.column_indices, dtype=np.int64)
        offsets = np.asarray(self.row_offsets, dtype=np.int64)
        if data.ndim != 1 or not np.all(np.isfinite(data)):
            raise ValueError("data must be a finite vector")
        if columns.ndim != 1 or offsets.ndim != 1:
            raise ValueError("CSR indices and offsets must be vectors")
        if not (
            isinstance(self.shape, tuple)
            and len(self.shape) == 2
            and all(type(dimension) is int for dimension in self.shape)
        ):
            raise TypeError("shape must be a pair of built-in integers")
        rows, columns_count = self.shape
        if rows < 0 or columns_count < 0:
            raise ValueError("shape dimensions must be nonnegative")
        if offsets.shape != (rows + 1,):
            raise ValueError("row_offsets length must equal rows plus one")
        if columns.shape != data.shape:
            raise ValueError("column_indices and data must have equal lengths")
        if np.any(data == 0.0):
            raise ValueError("canonical CSR data must not store explicit zeros")
        if offsets[0] != 0 or offsets[-1] != data.size:
            raise ValueError("row_offsets must span all stored entries")
        if np.any(offsets[1:] < offsets[:-1]):
            raise ValueError("row_offsets must be nondecreasing")
        if np.any(columns < 0) or np.any(columns >= columns_count):
            raise ValueError("column indices must lie within the matrix shape")
        for row in range(rows):
            row_columns = columns[offsets[row] : offsets[row + 1]]
            if np.any(row_columns[1:] <= row_columns[:-1]):
                raise ValueError(
                    "canonical CSR column indices must increase within each row"
                )
        if not isinstance(self.unit, PhysicalUnit | Unitless):
            raise TypeError("unit must be PhysicalUnit or Unitless")
        immutable_data = np.frombuffer(data.astype("<f8").tobytes(), dtype="<f8")
        immutable_columns = np.frombuffer(columns.astype("<i8").tobytes(), dtype="<i8")
        immutable_offsets = np.frombuffer(offsets.astype("<i8").tobytes(), dtype="<i8")
        object.__setattr__(self, "data", immutable_data)
        object.__setattr__(self, "column_indices", immutable_columns)
        object.__setattr__(self, "row_offsets", immutable_offsets)

    @classmethod
    def from_csr(
        cls,
        magnitude: SparseMatrix,
        unit: ModelSystemUnit,
    ) -> SparseMatrixQuantity:
        """Canonicalize one SciPy sparse matrix as an immutable CSR quantity."""
        if not isinstance(magnitude, sparse.spmatrix | sparse.sparray):
            raise TypeError("magnitude must be a SciPy sparse matrix or sparse array")
        if magnitude.dtype.kind not in "iuf":
            raise TypeError(
                "magnitude must contain real numeric values excluding booleans"
            )
        canonical = sparse.csr_array(magnitude, dtype=np.float64)
        canonical.sum_duplicates()
        canonical.eliminate_zeros()
        canonical.sort_indices()
        return cls(
            data=np.asarray(canonical.data, dtype=np.float64),
            column_indices=np.asarray(canonical.indices, dtype=np.int64),
            row_offsets=np.asarray(canonical.indptr, dtype=np.int64),
            shape=(int(canonical.shape[0]), int(canonical.shape[1])),
            unit=unit,
        )

    @property
    def nonzero_count(self) -> int:
        """Return the number of explicitly stored nonzero entries."""
        return int(self.data.size)

    def to_csr(self) -> sparse.csr_array:
        """Return a fresh SciPy CSR array with the represented values."""
        return sparse.csr_array(
            (
                self.data.copy(),
                self.column_indices.copy(),
                self.row_offsets.copy(),
            ),
            shape=self.shape,
            dtype=np.float64,
        )

    def to_dense(self) -> MatrixQuantity:
        """Materialize the represented matrix at an explicit dense boundary."""
        return MatrixQuantity(np.asarray(self.to_csr().toarray()), self.unit)


@dataclass(frozen=True, slots=True, eq=False)
class ComplexSparseMatrixQuantity(ModelSystemQuantity):
    """Represent an immutable finite complex128 matrix in canonical CSR storage.

    The quantity owns canonical CSR arrays backed by immutable bytes. Duplicate entries
    are summed, explicit zeros are eliminated, and indices are sorted by
    :meth:`from_csr`. Direct construction accepts only already-canonical CSR arrays.
    :meth:`to_csr` and :meth:`to_dense` are explicit fresh-copy boundaries, preventing
    mutation of maintained storage and preventing implicit densification.
    """

    data: ComplexVector
    column_indices: IntegerVector
    row_offsets: IntegerVector
    shape: tuple[int, int]
    unit: ModelSystemUnit

    def __post_init__(self) -> None:
        """Validate canonical complex CSR structure and install immutable storage."""
        if not isinstance(self.data, np.ndarray):
            raise TypeError("data must be a numpy.ndarray")
        if not isinstance(self.column_indices, np.ndarray):
            raise TypeError("column_indices must be a numpy.ndarray")
        if not isinstance(self.row_offsets, np.ndarray):
            raise TypeError("row_offsets must be a numpy.ndarray")
        if self.data.dtype.kind not in "iufc":
            raise TypeError("data must contain numeric values excluding booleans")
        if self.column_indices.dtype.kind not in "iu":
            raise TypeError("column_indices must contain integers excluding booleans")
        if self.row_offsets.dtype.kind not in "iu":
            raise TypeError("row_offsets must contain integers excluding booleans")
        data = np.asarray(self.data, dtype=np.complex128)
        columns = np.asarray(self.column_indices, dtype=np.int64)
        offsets = np.asarray(self.row_offsets, dtype=np.int64)
        if data.ndim != 1:
            raise ValueError("data must be a vector")
        if not np.all(np.isfinite(data.real)) or not np.all(np.isfinite(data.imag)):
            raise ValueError("data must contain finite complex values")
        if columns.ndim != 1 or offsets.ndim != 1:
            raise ValueError("CSR indices and offsets must be vectors")
        if not (
            isinstance(self.shape, tuple)
            and len(self.shape) == 2
            and all(type(dimension) is int for dimension in self.shape)
        ):
            raise TypeError("shape must be a pair of built-in integers")
        rows, columns_count = self.shape
        if rows < 0 or columns_count < 0:
            raise ValueError("shape dimensions must be nonnegative")
        if offsets.shape != (rows + 1,):
            raise ValueError("row_offsets length must equal rows plus one")
        if columns.shape != data.shape:
            raise ValueError("column_indices and data must have equal lengths")
        if np.any(data == 0.0 + 0.0j):
            raise ValueError("canonical CSR data must not store explicit zeros")
        if offsets[0] != 0 or offsets[-1] != data.size:
            raise ValueError("row_offsets must span all stored entries")
        if np.any(offsets[1:] < offsets[:-1]):
            raise ValueError("row_offsets must be nondecreasing")
        if np.any(columns < 0) or np.any(columns >= columns_count):
            raise ValueError("column indices must lie within the matrix shape")
        for row in range(rows):
            row_columns = columns[offsets[row] : offsets[row + 1]]
            if np.any(row_columns[1:] <= row_columns[:-1]):
                raise ValueError(
                    "canonical CSR column indices must increase within each row"
                )
        if not isinstance(self.unit, PhysicalUnit | Unitless):
            raise TypeError("unit must be PhysicalUnit or Unitless")
        immutable_data = np.frombuffer(data.astype("<c16").tobytes(), dtype="<c16")
        immutable_columns = np.frombuffer(columns.astype("<i8").tobytes(), dtype="<i8")
        immutable_offsets = np.frombuffer(offsets.astype("<i8").tobytes(), dtype="<i8")
        object.__setattr__(self, "data", immutable_data)
        object.__setattr__(self, "column_indices", immutable_columns)
        object.__setattr__(self, "row_offsets", immutable_offsets)

    @classmethod
    def from_csr(
        cls,
        magnitude: SparseMatrix,
        unit: ModelSystemUnit,
    ) -> ComplexSparseMatrixQuantity:
        """Canonicalize a numeric SciPy sparse matrix as complex128 CSR."""
        if not isinstance(magnitude, sparse.spmatrix | sparse.sparray):
            raise TypeError("magnitude must be a SciPy sparse matrix or sparse array")
        if magnitude.dtype.kind not in "iufc":
            raise TypeError("magnitude must contain numeric values excluding booleans")
        canonical = sparse.csr_array(magnitude, dtype=np.complex128)
        canonical.sum_duplicates()
        canonical.eliminate_zeros()
        canonical.sort_indices()
        return cls(
            data=np.asarray(canonical.data, dtype=np.complex128),
            column_indices=np.asarray(canonical.indices, dtype=np.int64),
            row_offsets=np.asarray(canonical.indptr, dtype=np.int64),
            shape=(int(canonical.shape[0]), int(canonical.shape[1])),
            unit=unit,
        )

    @property
    def nonzero_count(self) -> int:
        """Return the number of explicitly stored nonzero entries."""
        return int(self.data.size)

    def to_csr(self) -> sparse.csr_array:
        """Return a fresh SciPy complex128 CSR array."""
        return sparse.csr_array(
            (
                self.data.copy(),
                self.column_indices.copy(),
                self.row_offsets.copy(),
            ),
            shape=self.shape,
            dtype=np.complex128,
        )

    def to_dense(self) -> ComplexMatrixQuantity:
        """Materialize the represented matrix at an explicit dense boundary."""
        return ComplexMatrixQuantity(np.asarray(self.to_csr().toarray()), self.unit)


class PintUnitConverter:
    """Own Pint-backed validation, dimensional compatibility, and conversion."""

    __slots__ = ("registry",)

    registry: pint.UnitRegistry

    def __init__(self) -> None:
        """Use the one fixed registry that owns model-system conversions."""
        self.registry = MODEL_SYSTEM_PINT_REGISTRY

    def validate(self, unit: ModelSystemUnit) -> None:
        """Validate that one project unit is a parseable Pint unit."""
        if not isinstance(unit, PhysicalUnit | Unitless):
            raise TypeError("unit must be PhysicalUnit or Unitless")
        parsed = self.registry.Unit(unit.expression)
        if isinstance(unit, Unitless) and not parsed.dimensionless:
            raise ValueError("Unitless must be dimensionless")
        if isinstance(unit, PhysicalUnit) and parsed.dimensionless:
            raise ValueError("PhysicalUnit must not be dimensionless")

    def compatible(self, left: ModelSystemUnit, right: ModelSystemUnit) -> bool:
        """Return whether two validated units have equal dimensionality."""
        self.validate(left)
        self.validate(right)
        left_unit = self.registry.Unit(left.expression)
        right_unit = self.registry.Unit(right.expression)
        return bool(left_unit.is_compatible_with(right_unit))

    def conversion_factor(
        self, source: ModelSystemUnit, target: ModelSystemUnit
    ) -> float:
        """Return the multiplicative factor from ``source`` to ``target``.

        Raises
        ------
        ValueError
            If the units have incompatible dimensionality or the conversion has a
            nonzero affine offset. Scalar and real-vector conversion methods remain the
            explicit boundaries for affine conversions.
        """
        self.validate(source)
        self.validate(target)
        if not self.compatible(source, target):
            raise ValueError("source and target units are dimensionally incompatible")
        zero = self.registry.Quantity(0.0, source.expression).to(target.expression)
        if float(zero.magnitude) != 0.0:
            raise ValueError(
                "conversion_factor requires a conversion with zero affine offset"
            )
        quantity = self.registry.Quantity(1.0, source.expression)
        return float(quantity.to(target.expression).magnitude)

    def convert_scalar(
        self, quantity: ScalarQuantity, target: ModelSystemUnit
    ) -> ScalarQuantity:
        """Convert one scalar quantity, including affine unit conversions."""
        if not isinstance(quantity, ScalarQuantity):
            raise TypeError("quantity must be ScalarQuantity")
        self.validate(target)
        try:
            converted = self.registry.Quantity(
                quantity.magnitude,
                quantity.unit.expression,
            ).to(target.expression)
        except pint.errors.DimensionalityError as error:
            raise ValueError(
                "source and target units are dimensionally incompatible"
            ) from error
        return ScalarQuantity(
            magnitude=float(converted.magnitude),
            unit=target,
        )

    def convert_vector(
        self, quantity: VectorQuantity, target: ModelSystemUnit
    ) -> VectorQuantity:
        """Convert one vector quantity, including affine unit conversions."""
        if not isinstance(quantity, VectorQuantity):
            raise TypeError("quantity must be VectorQuantity")
        self.validate(target)
        try:
            converted = self.registry.Quantity(
                quantity.magnitude,
                quantity.unit.expression,
            ).to(target.expression)
        except pint.errors.DimensionalityError as error:
            raise ValueError(
                "source and target units are dimensionally incompatible"
            ) from error
        return VectorQuantity(
            magnitude=np.asarray(converted.magnitude, dtype=np.float64),
            unit=target,
        )

    def convert_matrix(
        self, quantity: MatrixQuantity, target: ModelSystemUnit
    ) -> MatrixQuantity:
        """Convert one dense real matrix through a multiplicative unit factor."""
        if not isinstance(quantity, MatrixQuantity):
            raise TypeError("quantity must be MatrixQuantity")
        factor = self.conversion_factor(quantity.unit, target)
        return MatrixQuantity(quantity.magnitude * factor, target)

    def convert_complex_vector(
        self, quantity: ComplexVectorQuantity, target: ModelSystemUnit
    ) -> ComplexVectorQuantity:
        """Convert one complex vector through a multiplicative unit factor."""
        if not isinstance(quantity, ComplexVectorQuantity):
            raise TypeError("quantity must be ComplexVectorQuantity")
        factor = self.conversion_factor(quantity.unit, target)
        return ComplexVectorQuantity(quantity.magnitude * factor, target)

    def convert_complex_matrix(
        self, quantity: ComplexMatrixQuantity, target: ModelSystemUnit
    ) -> ComplexMatrixQuantity:
        """Convert one dense complex matrix through a multiplicative unit factor."""
        if not isinstance(quantity, ComplexMatrixQuantity):
            raise TypeError("quantity must be ComplexMatrixQuantity")
        factor = self.conversion_factor(quantity.unit, target)
        return ComplexMatrixQuantity(quantity.magnitude * factor, target)

    def convert_sparse_matrix(
        self, quantity: SparseMatrixQuantity, target: ModelSystemUnit
    ) -> SparseMatrixQuantity:
        """Convert one sparse real matrix through a multiplicative unit factor."""
        if not isinstance(quantity, SparseMatrixQuantity):
            raise TypeError("quantity must be SparseMatrixQuantity")
        factor = self.conversion_factor(quantity.unit, target)
        return SparseMatrixQuantity.from_csr(quantity.to_csr() * factor, target)

    def convert_complex_sparse_matrix(
        self, quantity: ComplexSparseMatrixQuantity, target: ModelSystemUnit
    ) -> ComplexSparseMatrixQuantity:
        """Convert one sparse complex matrix through a multiplicative unit factor."""
        if not isinstance(quantity, ComplexSparseMatrixQuantity):
            raise TypeError("quantity must be ComplexSparseMatrixQuantity")
        factor = self.conversion_factor(quantity.unit, target)
        return ComplexSparseMatrixQuantity.from_csr(quantity.to_csr() * factor, target)


ModelSystemUnitConverter = PintUnitConverter()
"""Shared converter using one fixed Pint registry for model-system operations."""
