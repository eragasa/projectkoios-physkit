"""Range-family aggregate result invariant tests."""

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
    BlockHoppingRangeBandEvaluation,
    BlockHoppingRangeBandEvaluator,
    BlockHoppingRangeBandSelection,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    ScalarQuantity,
    Unitless,
)


def _model() -> BlockHoppingModel:
    dimension = LatticeDimension.ONE
    return BlockHoppingModel(
        "source",
        dimension,
        tuple(
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (value,)),
                ComplexMatrixQuantity(np.array([[coefficient]]), Unitless()),
            )
            for value, coefficient in ((-1, 0.5), (0, 2.0), (1, 0.5))
        ),
        Unitless(),
        "zero",
        "basis",
    )


def _selection(maximum: int, identifier: str) -> BlockHoppingRangeBandSelection:
    return BlockHoppingRangeBandSelection(
        tuple(
            LatticeDisplacement(LatticeDimension.ONE, (value,))
            for value in range(-maximum, maximum + 1)
        ),
        identifier,
    )


def test__init__retains_source_entries_and_derived_inputs() -> None:
    source = _model()
    coordinates = MatrixQuantity(np.array([[0.0], [0.25]]), Unitless())
    tolerance = BlockHoppingBandTolerance(ScalarQuantity(0.0, Unitless()), 0.0)
    result = BlockHoppingRangeBandEvaluator().evaluate(
        source,
        (_selection(0, "onsite"), _selection(1, "full retained")),
        coordinates,
        tolerance,
    )

    rebuilt = BlockHoppingRangeBandEvaluation(
        result.source_band_evaluation,
        result.entries,
    )

    assert rebuilt.source is source
    assert rebuilt.fractional_wave_vectors is coordinates
    assert rebuilt.tolerance is tolerance
    assert rebuilt.entries is result.entries


def test__init__rejects_wrong_empty_or_duplicate_entries() -> None:
    source = _model()
    coordinates = MatrixQuantity(np.array([[0.0]]), Unitless())
    tolerance = BlockHoppingBandTolerance(ScalarQuantity(0.0, Unitless()), 0.0)
    result = BlockHoppingRangeBandEvaluator().evaluate(
        source,
        (_selection(0, "onsite"), _selection(1, "full retained")),
        coordinates,
        tolerance,
    )

    with pytest.raises(TypeError, match="source_band_evaluation"):
        BlockHoppingRangeBandEvaluation(cast(Any, object()), result.entries)
    with pytest.raises(TypeError, match="entries must be a tuple"):
        BlockHoppingRangeBandEvaluation(
            result.source_band_evaluation,
            cast(Any, list(result.entries)),
        )
    with pytest.raises(ValueError, match="entries must be nonempty"):
        BlockHoppingRangeBandEvaluation(result.source_band_evaluation, ())
    with pytest.raises(ValueError, match="identifiers must be unique"):
        BlockHoppingRangeBandEvaluation(
            result.source_band_evaluation,
            (result.entries[0], result.entries[0]),
        )


def test__init__rejects_duplicate_inventory_or_unshared_coordinates() -> None:
    source = _model()
    coordinates = MatrixQuantity(np.array([[0.0]]), Unitless())
    tolerance = BlockHoppingBandTolerance(ScalarQuantity(0.0, Unitless()), 0.0)
    evaluator = BlockHoppingRangeBandEvaluator()
    first = evaluator.evaluate(
        source,
        (_selection(0, "onsite-a"),),
        coordinates,
        tolerance,
    )
    same_inventory = evaluator.evaluate(
        source,
        (_selection(0, "onsite-b"),),
        coordinates,
        tolerance,
    )
    different_coordinates = evaluator.evaluate(
        source,
        (_selection(1, "full"),),
        MatrixQuantity(np.array([[0.0]]), Unitless()),
        tolerance,
    )

    with pytest.raises(ValueError, match="inventories must be unique"):
        BlockHoppingRangeBandEvaluation(
            first.source_band_evaluation,
            (first.entries[0], same_inventory.entries[0]),
        )
    with pytest.raises(ValueError, match="shared coordinates"):
        BlockHoppingRangeBandEvaluation(
            first.source_band_evaluation,
            (different_coordinates.entries[0],),
        )
