"""Correlation tests for raw block-power results."""

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
from projectkoios.physkit.solidstate.tb.diagnostics import (
    BlockHoppingPower,
    BlockHoppingPowerAnalyzer,
    BlockHoppingPowerResult,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
)


def _model() -> BlockHoppingModel:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    return BlockHoppingModel(
        "model",
        dimension,
        (
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (0,)),
                ComplexMatrixQuantity(np.array([[2.0]]), unit),
            ),
        ),
        unit,
        "zero",
        "basis",
    )


def test__init__accepts_only_canonical_power_inventory() -> None:
    canonical = BlockHoppingPowerAnalyzer().execute(_model())

    rebuilt = BlockHoppingPowerResult(canonical.model, canonical.powers)

    assert rebuilt.powers == canonical.powers


def test__init__rejects_wrong_types_and_noncanonical_values() -> None:
    model = _model()
    canonical = BlockHoppingPowerAnalyzer().execute(model)

    with pytest.raises(TypeError, match="model must be BlockHoppingModel"):
        BlockHoppingPowerResult(cast(Any, object()), canonical.powers)
    with pytest.raises(TypeError, match="BlockHoppingPower"):
        BlockHoppingPowerResult(model, cast(Any, (object(),)))
    wrong = BlockHoppingPower(
        canonical.powers[0].displacement,
        ScalarQuantity(5.0, PhysicalUnit("(electron_volt) ** 2")),
    )
    with pytest.raises(ValueError, match="canonical Frobenius powers"):
        BlockHoppingPowerResult(model, (wrong,))
