"""Construction tests for matrix-valued hopping models."""

from dataclasses import FrozenInstanceError
from typing import Any, cast

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    BlockHoppingModel,
    BlockHoppingTerm,
)
from projectkoios.physkit.units import ComplexMatrixQuantity, PhysicalUnit, Unitless


def _term(
    dimension: LatticeDimension,
    components: tuple[int, ...],
    values: np.ndarray,
    unit: PhysicalUnit | Unitless,
) -> BlockHoppingTerm:
    closed_components: tuple[int] | tuple[int, int] | tuple[int, int, int]
    if dimension is LatticeDimension.ONE:
        closed_components = (components[0],)
    elif dimension is LatticeDimension.TWO:
        closed_components = (components[0], components[1])
    else:
        closed_components = (components[0], components[1], components[2])
    return BlockHoppingTerm(
        LatticeDisplacement(dimension, closed_components),
        ComplexMatrixQuantity(values, unit),
    )


def test__init__retains_ordered_basis_and_energy_identity() -> None:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    negative = _term(dimension, (-1,), np.eye(2), unit)
    onsite = _term(dimension, (0,), np.diag([2.0, 3.0]), unit)
    positive = _term(dimension, (1,), np.eye(2), unit)

    model = BlockHoppingModel(
        identifier="two-band chain",
        dimension=dimension,
        terms=(negative, onsite, positive),
        energy_unit=unit,
        energy_reference="common zero",
        basis_identifier="ordered basis: a, b",
    )

    assert model.terms == (negative, onsite, positive)
    assert model.energy_unit is unit
    assert model.energy_reference == "common zero"
    assert model.basis_identifier == "ordered basis: a, b"
    assert model.basis_size == 2
    with pytest.raises(FrozenInstanceError):
        model.identifier = "other"  # type: ignore[misc]


def test__init__accepts_explicit_unitless_blocks() -> None:
    dimension = LatticeDimension.ONE
    unit = Unitless()
    term = _term(dimension, (0,), np.ones((1, 1)), unit)

    model = BlockHoppingModel(
        "dimensionless",
        dimension,
        (term,),
        unit,
        "zero",
        "one basis state",
    )

    assert model.basis_size == 1
    assert isinstance(model.energy_unit, Unitless)


def test__init__rejects_inconsistent_dimension_shape_and_unit() -> None:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    valid = _term(dimension, (0,), np.eye(2), unit)
    wrong_dimension = _term(LatticeDimension.TWO, (1, 0), np.eye(2), unit)
    wrong_shape = _term(dimension, (1,), np.eye(3), unit)
    wrong_unit = _term(dimension, (1,), np.eye(2), PhysicalUnit("joule"))

    for term, message in (
        (wrong_dimension, "dimensions must match"),
        (wrong_shape, "common square shape"),
        (wrong_unit, "units must match"),
    ):
        with pytest.raises(ValueError, match=message):
            BlockHoppingModel(
                "invalid",
                dimension,
                (valid, term),
                unit,
                "zero",
                "basis",
            )


def test__init__rejects_duplicate_or_unsorted_displacements() -> None:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    negative = _term(dimension, (-1,), np.eye(2), unit)
    positive = _term(dimension, (1,), np.eye(2), unit)

    for terms in ((negative, negative), (positive, negative)):
        with pytest.raises(ValueError, match="sorted and unique"):
            BlockHoppingModel(
                "invalid",
                dimension,
                terms,
                unit,
                "zero",
                "basis",
            )


def test__init__rejects_empty_or_nonowned_contract_values() -> None:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    term = _term(dimension, (0,), np.eye(2), unit)

    with pytest.raises(ValueError, match="terms must be nonempty"):
        BlockHoppingModel("empty", dimension, (), unit, "zero", "basis")
    with pytest.raises(TypeError, match="terms must be a tuple"):
        BlockHoppingModel(
            "list",
            dimension,
            cast(Any, [term]),
            unit,
            "zero",
            "basis",
        )
    with pytest.raises(TypeError, match="BlockHoppingTerm"):
        BlockHoppingModel(
            "wrong term",
            dimension,
            (cast(Any, object()),),
            unit,
            "zero",
            "basis",
        )
    with pytest.raises(TypeError, match="PhysicalUnit or Unitless"):
        BlockHoppingModel(
            "wrong unit",
            dimension,
            (term,),
            cast(Any, "electron_volt"),
            "zero",
            "basis",
        )
    with pytest.raises(ValueError, match="basis_identifier must be nonempty"):
        BlockHoppingModel("empty basis", dimension, (term,), unit, "zero", "")
