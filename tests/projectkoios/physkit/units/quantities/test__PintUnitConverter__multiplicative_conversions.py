"""Verification of multiplicative dense and sparse quantity conversions."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.units.quantities import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    ComplexVectorQuantity,
    MatrixQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    SparseMatrixQuantity,
)


def test__conversion_factor__returns_scale_for_multiplicative_units() -> None:
    """A factor is available only where zero maps to zero."""
    assert (
        ModelSystemUnitConverter.conversion_factor(
            PhysicalUnit("meter"),
            PhysicalUnit("centimeter"),
        )
        == 100.0
    )


def test__conversion_factor__rejects_affine_and_incompatible_units() -> None:
    """A single factor must not silently erase offsets or dimensional mismatch."""
    with pytest.raises(ValueError, match="zero affine offset"):
        ModelSystemUnitConverter.conversion_factor(
            PhysicalUnit("degree_Celsius"),
            PhysicalUnit("kelvin"),
        )
    with pytest.raises(ValueError, match="dimensionally incompatible"):
        ModelSystemUnitConverter.conversion_factor(
            PhysicalUnit("meter"),
            PhysicalUnit("second"),
        )


def test__dense_converters__preserve_real_and_complex_values() -> None:
    """Dense linear-algebra quantities use the declared multiplicative factor."""
    meter = PhysicalUnit("meter")
    centimeter = PhysicalUnit("centimeter")

    real = ModelSystemUnitConverter.convert_matrix(
        MatrixQuantity(np.asarray([[1.0, -2.0]]), meter),
        centimeter,
    )
    vector = ModelSystemUnitConverter.convert_complex_vector(
        ComplexVectorQuantity(np.asarray([1.0 + 2.0j]), meter),
        centimeter,
    )
    matrix = ModelSystemUnitConverter.convert_complex_matrix(
        ComplexMatrixQuantity(np.asarray([[1.0 - 2.0j]]), meter),
        centimeter,
    )

    np.testing.assert_array_equal(real.magnitude, [[100.0, -200.0]])
    np.testing.assert_array_equal(vector.magnitude, [100.0 + 200.0j])
    np.testing.assert_array_equal(matrix.magnitude, [[100.0 - 200.0j]])
    assert real.unit == vector.unit == matrix.unit == centimeter


def test__sparse_converters__preserve_canonical_storage_without_densifying() -> None:
    """Sparse conversion scales stored values and keeps canonical CSR results."""
    meter = PhysicalUnit("meter")
    millimeter = PhysicalUnit("millimeter")
    real_source = SparseMatrixQuantity.from_csr(
        sparse.csr_array([[0.0, 2.0], [-3.0, 0.0]]),
        meter,
    )
    complex_source = ComplexSparseMatrixQuantity.from_csr(
        sparse.csr_array([[0.0j, 1.0 + 2.0j]]),
        meter,
    )

    real = ModelSystemUnitConverter.convert_sparse_matrix(
        real_source,
        millimeter,
    )
    complex_result = ModelSystemUnitConverter.convert_complex_sparse_matrix(
        complex_source,
        millimeter,
    )

    np.testing.assert_array_equal(
        real.to_csr().toarray(), [[0.0, 2000.0], [-3000.0, 0.0]]
    )
    np.testing.assert_array_equal(
        complex_result.to_csr().toarray(),
        [[0.0j, 1000.0 + 2000.0j]],
    )
    assert real.nonzero_count == real_source.nonzero_count
    assert complex_result.nonzero_count == complex_source.nonzero_count
    assert real.unit == complex_result.unit == millimeter


@pytest.mark.parametrize(
    "convert, quantity",
    (
        (
            ModelSystemUnitConverter.convert_matrix,
            MatrixQuantity(
                np.asarray([[0.0]]),
                PhysicalUnit("degree_Celsius"),
            ),
        ),
        (
            ModelSystemUnitConverter.convert_complex_vector,
            ComplexVectorQuantity(
                np.asarray([0.0j]),
                PhysicalUnit("degree_Celsius"),
            ),
        ),
        (
            ModelSystemUnitConverter.convert_complex_matrix,
            ComplexMatrixQuantity(
                np.asarray([[0.0j]]),
                PhysicalUnit("degree_Celsius"),
            ),
        ),
        (
            ModelSystemUnitConverter.convert_sparse_matrix,
            SparseMatrixQuantity.from_csr(
                sparse.csr_array([[1.0]]),
                PhysicalUnit("degree_Celsius"),
            ),
        ),
        (
            ModelSystemUnitConverter.convert_complex_sparse_matrix,
            ComplexSparseMatrixQuantity.from_csr(
                sparse.csr_array([[1.0j]]),
                PhysicalUnit("degree_Celsius"),
            ),
        ),
    ),
)
def test__linear_algebra_converters__reject_affine_units(convert, quantity) -> None:
    """Offset conversion cannot be represented by one linear scale factor."""
    with pytest.raises(ValueError, match="zero affine offset"):
        convert(quantity, PhysicalUnit("kelvin"))
