"""Tests for correlated scalar-hopping Hermiticity results."""

import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    ScalarHoppingModel,
    ScalarHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.diagnostics import (
    ScalarHoppingCoefficientComparison,
    ScalarHoppingHermiticityResult,
)
from projectkoios.physkit.units import Unitless


def test__init__rejects_comparisons_that_do_not_match_the_model() -> None:
    dimension = LatticeDimension.ONE
    source = LatticeDisplacement(dimension, (1,))
    target = LatticeDisplacement(dimension, (-1,))
    model = ScalarHoppingModel(
        "one_sided",
        dimension,
        (ScalarHoppingTerm(source, 1.0, 0.0),),
        Unitless(),
        "zero",
        "basis",
    )
    incorrect = ScalarHoppingCoefficientComparison(source, target, 0.0)

    with pytest.raises(ValueError, match="canonical Hermiticity comparisons"):
        ScalarHoppingHermiticityResult(model, (incorrect,))
