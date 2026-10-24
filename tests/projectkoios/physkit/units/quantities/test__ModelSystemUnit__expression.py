"""Static and runtime ownership tests for nominal unit expressions."""

from typing import assert_type, get_type_hints

import pytest

from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    ComplexVectorQuantity,
    MatrixQuantity,
    ModelSystemUnit,
    PhysicalUnit,
    ScalarQuantity,
    SparseMatrixQuantity,
    Unitless,
    VectorQuantity,
)


def _nominal_expression(unit: ModelSystemUnit) -> str:
    """Exercise the expression contract through its nominal static type."""
    return assert_type(unit.expression, str)


@pytest.mark.parametrize(
    ("unit", "expected"),
    [
        (PhysicalUnit("electron_volt"), "electron_volt"),
        (Unitless(), "dimensionless"),
    ],
)
def test__expression__is_shared_by_concrete_nominal_units(
    unit: ModelSystemUnit,
    expected: str,
) -> None:
    assert isinstance(unit, ModelSystemUnit)
    assert _nominal_expression(unit) == expected


def test__expression__adds_no_runtime_descriptor_to_nominal_base() -> None:
    assert "expression" not in ModelSystemUnit.__dict__


@pytest.mark.parametrize(
    "quantity_type",
    [
        ScalarQuantity,
        VectorQuantity,
        ComplexVectorQuantity,
        MatrixQuantity,
        ComplexMatrixQuantity,
        SparseMatrixQuantity,
        ComplexSparseMatrixQuantity,
    ],
)
def test__unit_field__retains_canonical_nominal_annotation(
    quantity_type: type[object],
) -> None:
    assert get_type_hints(quantity_type)["unit"] is ModelSystemUnit
