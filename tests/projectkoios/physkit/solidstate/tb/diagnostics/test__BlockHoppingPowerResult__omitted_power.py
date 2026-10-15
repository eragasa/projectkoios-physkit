"""Explicit omitted-power inventory tests."""

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
from projectkoios.physkit.solidstate.tb.diagnostics import BlockHoppingPowerAnalyzer
from projectkoios.physkit.units import ComplexMatrixQuantity, PhysicalUnit


def _result():
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    model = BlockHoppingModel(
        "power model",
        dimension,
        tuple(
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (representative,)),
                ComplexMatrixQuantity(np.array([[value]]), unit),
            )
            for representative, value in ((-1, 2.0), (0, 3.0), (1, 4.0))
        ),
        unit,
        "zero",
        "basis",
    )
    return BlockHoppingPowerAnalyzer().execute(model)


def test__omitted_power__sums_outside_explicit_inventory() -> None:
    result = _result()

    omitted = result.omitted_power((result.model.terms[1].displacement,))
    all_omitted = result.omitted_power(())
    none_omitted = result.omitted_power(
        tuple(term.displacement for term in result.model.terms)
    )

    assert omitted.magnitude == 20.0
    assert all_omitted.magnitude == 29.0
    assert none_omitted.magnitude == 0.0
    assert omitted.unit == PhysicalUnit("(electron_volt) ** 2")


def test__omitted_power__rejects_invalid_inventories() -> None:
    result = _result()
    negative = result.model.terms[0].displacement
    zero = result.model.terms[1].displacement

    with pytest.raises(TypeError, match="tuple of LatticeDisplacement"):
        result.omitted_power(cast(Any, [zero]))
    with pytest.raises(ValueError, match="sorted and unique"):
        result.omitted_power((zero, negative))
    with pytest.raises(ValueError, match="sorted and unique"):
        result.omitted_power((zero, zero))
    with pytest.raises(ValueError, match="belong to the model"):
        result.omitted_power((LatticeDisplacement(LatticeDimension.ONE, (2,)),))
    with pytest.raises(ValueError, match="belong to the model"):
        result.omitted_power((LatticeDisplacement(LatticeDimension.TWO, (0, 0)),))
