"""Tests for scalar-hopping comparisons under supplied lattice operations."""

import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    ScalarHoppingModel,
    ScalarHoppingTerm,
)
from projectkoios.physkit.periodic.lattice.symmetry import IntegralLatticeOperation
from projectkoios.physkit.solidstate.tb.diagnostics import (
    ScalarHoppingSymmetryAnalyzer,
)
from projectkoios.physkit.units import PhysicalUnit


def _two_dimensional_model() -> ScalarHoppingModel:
    dimension = LatticeDimension.TWO
    return ScalarHoppingModel(
        "symmetry_diagnostic",
        dimension,
        (
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (-1, 0)),
                2.0,
                0.0,
            ),
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (0, -1)),
                3.0,
                0.0,
            ),
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (0, 0)),
                4.0,
                1.0,
            ),
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (0, 1)),
                1.0,
                0.0,
            ),
        ),
        PhysicalUnit("electron_volt"),
        "zero",
        "basis",
    )


def test__execute__compares_every_source_with_its_transformed_target() -> None:
    model = _two_dimensional_model()
    swap = IntegralLatticeOperation(
        "swap",
        LatticeDimension.TWO,
        ((0, 1), (1, 0)),
    )

    result = ScalarHoppingSymmetryAnalyzer().execute(model, swap)

    assert result.model is model
    assert result.operation is swap
    assert tuple(
        (
            comparison.source.components,
            comparison.target.components,
            comparison.absolute_coefficient_difference,
        )
        for comparison in result.comparisons
    ) == (
        ((-1, 0), (0, -1), 1.0),
        ((0, -1), (-1, 0), 1.0),
        ((0, 0), (0, 0), 0.0),
        ((0, 1), (1, 0), None),
    )
    assert tuple(target.components for target in result.missing_targets) == ((1, 0),)
    assert result.maximum_absolute_coefficient_difference == 1.0


def test__execute__accepts_a_general_unimodular_operation() -> None:
    dimension = LatticeDimension.TWO
    model = ScalarHoppingModel(
        "shear_invariant_support",
        dimension,
        (
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (1, 0)),
                2.0,
                0.0,
            ),
        ),
        PhysicalUnit("electron_volt"),
        "zero",
        "basis",
    )
    shear = IntegralLatticeOperation(
        "shear",
        dimension,
        ((1, 1), (0, 1)),
    )

    result = ScalarHoppingSymmetryAnalyzer().execute(model, shear)

    assert result.maximum_absolute_coefficient_difference == 0.0
    assert not result.missing_targets


def test__execute__supports_three_dimensional_operations() -> None:
    dimension = LatticeDimension.THREE
    model = ScalarHoppingModel(
        "three_dimensional",
        dimension,
        (
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (0, 0, 1)),
                2.0,
                0.0,
            ),
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (0, 1, 0)),
                2.0,
                0.0,
            ),
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (1, 0, 0)),
                2.0,
                0.0,
            ),
        ),
        PhysicalUnit("electron_volt"),
        "zero",
        "basis",
    )
    axis_cycle = IntegralLatticeOperation(
        "axis_cycle",
        dimension,
        ((0, 0, 1), (1, 0, 0), (0, 1, 0)),
    )

    result = ScalarHoppingSymmetryAnalyzer().execute(model, axis_cycle)

    assert result.maximum_absolute_coefficient_difference == 0.0
    assert not result.missing_targets


def test__execute__rejects_a_dimensionally_incompatible_operation() -> None:
    model = _two_dimensional_model()
    one_dimensional_identity = IntegralLatticeOperation(
        "identity",
        LatticeDimension.ONE,
        ((1,),),
    )

    with pytest.raises(ValueError, match="dimension must match"):
        ScalarHoppingSymmetryAnalyzer().execute(
            model,
            one_dimensional_identity,
        )
