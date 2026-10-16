"""Range-family entry invariant tests."""

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
from projectkoios.physkit.solidstate.tb.bands import (
    BlockHoppingBandTolerance,
    BlockHoppingRangeBandEntry,
    BlockHoppingRangeBandEvaluator,
    BlockHoppingRangeBandSelection,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    ScalarQuantity,
    Unitless,
)


def _evaluation():
    dimension = LatticeDimension.ONE
    zero = LatticeDisplacement(dimension, (0,))
    source = BlockHoppingModel(
        "source",
        dimension,
        (BlockHoppingTerm(zero, ComplexMatrixQuantity(np.array([[2.0]]), Unitless())),),
        Unitless(),
        "zero",
        "basis",
    )
    return BlockHoppingRangeBandEvaluator().evaluate(
        source,
        (BlockHoppingRangeBandSelection((zero,), "retained"),),
        MatrixQuantity(np.array([[0.0]]), Unitless()),
        BlockHoppingBandTolerance(ScalarQuantity(0.0, Unitless()), 0.0),
    )


def test__init__retains_exact_truncation_and_band_evaluation() -> None:
    result = _evaluation()
    entry = result.entries[0]

    rebuilt = BlockHoppingRangeBandEntry(entry.truncation, entry.band_evaluation)

    assert rebuilt.truncation is entry.truncation
    assert rebuilt.band_evaluation is entry.band_evaluation


def test__init__rejects_wrong_types_or_model_identity() -> None:
    result = _evaluation()
    entry = result.entries[0]

    with pytest.raises(TypeError, match="truncation"):
        BlockHoppingRangeBandEntry(cast(Any, object()), entry.band_evaluation)
    with pytest.raises(TypeError, match="band_evaluation"):
        BlockHoppingRangeBandEntry(entry.truncation, cast(Any, object()))
    with pytest.raises(ValueError, match="exact truncated model"):
        BlockHoppingRangeBandEntry(
            entry.truncation,
            result.source_band_evaluation,
        )
