"""Construction tests for the block-hopping Bloch adapter."""

from typing import Any, cast

import pytest

from projectkoios.numerics.finite_periodic.finite_series import (
    FiniteFourierSeriesEvaluationDefinition,
    FiniteFourierSeriesEvaluationResult,
    FiniteFourierSeriesEvaluator,
)
from projectkoios.physkit.solidstate.tb.bloch import BlockHoppingBlochEvaluator


class _Evaluator(FiniteFourierSeriesEvaluator):
    """Minimal nominal evaluator for constructor verification."""

    def _evaluate_validated(
        self,
        definition: FiniteFourierSeriesEvaluationDefinition,
        /,
    ) -> FiniteFourierSeriesEvaluationResult:
        return FiniteFourierSeriesEvaluationResult(
            definition.dimension,
            definition.point_count,
            definition.term_count,
            (0.0j,) * definition.point_count,
        )


def test__init__retains_injected_nominal_evaluator() -> None:
    owner = _Evaluator()

    evaluator = BlockHoppingBlochEvaluator(owner)

    assert evaluator.finite_series_evaluator is owner


def test__init__constructs_default_nominal_evaluator() -> None:
    evaluator = BlockHoppingBlochEvaluator()

    assert isinstance(evaluator.finite_series_evaluator, FiniteFourierSeriesEvaluator)


def test__init__rejects_nonowned_evaluator() -> None:
    with pytest.raises(TypeError, match="FiniteFourierSeriesEvaluator or None"):
        BlockHoppingBlochEvaluator(cast(Any, object()))
