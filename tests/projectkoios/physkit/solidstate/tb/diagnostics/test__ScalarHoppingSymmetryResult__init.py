"""Tests for correlated scalar-hopping symmetry results."""

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
    ScalarHoppingCoefficientComparison,
    ScalarHoppingSymmetryResult,
)
from projectkoios.physkit.units import Unitless


def test__init__rejects_comparisons_that_do_not_match_the_operation() -> None:
    dimension = LatticeDimension.ONE
    displacement = LatticeDisplacement(dimension, (1,))
    model = ScalarHoppingModel(
        "one_term",
        dimension,
        (ScalarHoppingTerm(displacement, 1.0, 0.0),),
        Unitless(),
        "zero",
        "basis",
    )
    inversion = IntegralLatticeOperation("inversion", dimension, ((-1,),))
    incorrect = ScalarHoppingCoefficientComparison(
        displacement,
        displacement,
        0.0,
    )

    with pytest.raises(ValueError, match="supplied lattice operation"):
        ScalarHoppingSymmetryResult(model, inversion, (incorrect,))
