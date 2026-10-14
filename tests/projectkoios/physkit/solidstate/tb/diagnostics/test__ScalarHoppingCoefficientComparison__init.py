"""Tests for one scalar-hopping coefficient comparison."""

import math

import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.solidstate.tb.diagnostics import (
    ScalarHoppingCoefficientComparison,
)


def test__init__distinguishes_an_available_zero_difference_from_a_missing_target() -> (
    None
):
    source = LatticeDisplacement(LatticeDimension.TWO, (1, 0))
    target = LatticeDisplacement(LatticeDimension.TWO, (0, 1))

    available = ScalarHoppingCoefficientComparison(source, target, 0.0)
    missing = ScalarHoppingCoefficientComparison(source, target, None)

    assert available.target_is_present
    assert available.absolute_coefficient_difference == 0.0
    assert not missing.target_is_present
    assert missing.absolute_coefficient_difference is None


def test__init__rejects_incompatible_displacements_and_invalid_differences() -> None:
    one_dimensional = LatticeDisplacement(LatticeDimension.ONE, (0,))
    two_dimensional = LatticeDisplacement(LatticeDimension.TWO, (0, 0))

    with pytest.raises(ValueError, match="dimensions must agree"):
        ScalarHoppingCoefficientComparison(
            one_dimensional,
            two_dimensional,
            0.0,
        )
    with pytest.raises(TypeError, match="built-in float or None"):
        ScalarHoppingCoefficientComparison(
            one_dimensional,
            one_dimensional,
            0,  # type: ignore[arg-type]
        )
    with pytest.raises(ValueError, match="finite and nonnegative"):
        ScalarHoppingCoefficientComparison(
            one_dimensional,
            one_dimensional,
            math.inf,
        )
