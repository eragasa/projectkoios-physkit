"""Tests for scalar-hopping Hermiticity comparisons."""

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    ScalarHoppingModel,
    ScalarHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.diagnostics import (
    ScalarHoppingHermiticityAnalyzer,
)
from projectkoios.physkit.units import PhysicalUnit


def test__execute__compares_each_pair_once_and_identifies_missing_targets() -> None:
    dimension = LatticeDimension.TWO
    model = ScalarHoppingModel(
        "diagnostic_model",
        dimension,
        (
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (-1, 0)),
                2.0,
                -1.0,
            ),
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (0, 0)),
                3.0,
                0.25,
            ),
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (0, 1)),
                4.0,
                0.0,
            ),
            ScalarHoppingTerm(
                LatticeDisplacement(dimension, (1, 0)),
                2.0,
                1.5,
            ),
        ),
        PhysicalUnit("electron_volt"),
        "zero",
        "basis",
    )

    result = ScalarHoppingHermiticityAnalyzer().execute(model)

    assert result.model is model
    assert tuple(
        (
            comparison.source.components,
            comparison.target.components,
            comparison.absolute_coefficient_difference,
        )
        for comparison in result.comparisons
    ) == (
        ((-1, 0), (1, 0), 0.5),
        ((0, 0), (0, 0), 0.5),
        ((0, 1), (0, -1), None),
    )
    assert tuple(target.components for target in result.missing_targets) == ((0, -1),)
    assert result.maximum_absolute_coefficient_difference == 0.5


def test__execute__reports_no_maximum_when_every_target_is_missing() -> None:
    dimension = LatticeDimension.ONE
    model = ScalarHoppingModel(
        "one_sided",
        dimension,
        (ScalarHoppingTerm(LatticeDisplacement(dimension, (1,)), 1.0, 0.0),),
        PhysicalUnit("electron_volt"),
        "zero",
        "basis",
    )

    result = ScalarHoppingHermiticityAnalyzer().execute(model)

    assert result.maximum_absolute_coefficient_difference is None
    assert result.missing_targets[0].components == (-1,)


def test__execute__retains_an_exact_zero_difference() -> None:
    dimension = LatticeDimension.ONE
    model = ScalarHoppingModel(
        "hermitian",
        dimension,
        (
            ScalarHoppingTerm(LatticeDisplacement(dimension, (-1,)), 1.0, -2.0),
            ScalarHoppingTerm(LatticeDisplacement(dimension, (1,)), 1.0, 2.0),
        ),
        PhysicalUnit("electron_volt"),
        "zero",
        "basis",
    )

    result = ScalarHoppingHermiticityAnalyzer().execute(model)

    assert result.comparisons[0].target_is_present
    assert result.comparisons[0].absolute_coefficient_difference == 0.0
    assert result.maximum_absolute_coefficient_difference == 0.0
