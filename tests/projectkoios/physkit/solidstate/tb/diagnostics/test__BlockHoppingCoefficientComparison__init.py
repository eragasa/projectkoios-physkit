"""Invariant tests for matrix-hopping coefficient comparisons."""

from typing import Any, cast

import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.solidstate.tb.diagnostics import (
    BlockHoppingCoefficientComparison,
)
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity


def test__init__retains_present_and_missing_target_dispositions() -> None:
    source = LatticeDisplacement(LatticeDimension.ONE, (-1,))
    target = LatticeDisplacement(LatticeDimension.ONE, (1,))
    defect = ScalarQuantity(0.25, PhysicalUnit("electron_volt"))

    present = BlockHoppingCoefficientComparison(source, target, defect)
    missing = BlockHoppingCoefficientComparison(source, target, None)

    assert present.frobenius_defect is defect
    assert present.target_is_present
    assert missing.frobenius_defect is None
    assert not missing.target_is_present


def test__init__rejects_wrong_types_dimensions_and_negative_defect() -> None:
    one = LatticeDisplacement(LatticeDimension.ONE, (1,))
    two = LatticeDisplacement(LatticeDimension.TWO, (1, 0))

    with pytest.raises(TypeError, match="source must be LatticeDisplacement"):
        BlockHoppingCoefficientComparison(cast(Any, object()), one, None)
    with pytest.raises(TypeError, match="target must be LatticeDisplacement"):
        BlockHoppingCoefficientComparison(one, cast(Any, object()), None)
    with pytest.raises(ValueError, match="dimensions must agree"):
        BlockHoppingCoefficientComparison(one, two, None)
    with pytest.raises(TypeError, match="ScalarQuantity or None"):
        BlockHoppingCoefficientComparison(one, one, cast(Any, 0.0))
    with pytest.raises(ValueError, match="must be nonnegative"):
        BlockHoppingCoefficientComparison(
            one,
            one,
            ScalarQuantity(-1.0, PhysicalUnit("electron_volt")),
        )
