"""Construction tests for the block-hopping Bloch adapter."""

from typing import Any, cast

import pytest

from projectkoios.numerics.finite_periodic.batched_finite_series import (
    BatchedFiniteFourierSeriesEvaluationDefinition,
    BatchedFiniteFourierSeriesEvaluationResult,
    BatchedFiniteFourierSeriesEvaluator,
)
from projectkoios.numerics.finite_periodic.finite_series import (
    FiniteFourierSeriesEvaluationDefinition,
    FiniteFourierSeriesEvaluationResult,
    FiniteFourierSeriesEvaluator,
)
from projectkoios.physkit.solidstate.tb.bloch import BlockHoppingBlochEvaluator


class _BatchedEvaluator(BatchedFiniteFourierSeriesEvaluator):
    """Minimal shared-grid evaluator for constructor verification."""

    def _evaluate_validated(
        self,
        definition: BatchedFiniteFourierSeriesEvaluationDefinition,
        /,
    ) -> BatchedFiniteFourierSeriesEvaluationResult:
        return BatchedFiniteFourierSeriesEvaluationResult(
            definition.dimension,
            definition.point_count,
            definition.term_count,
            definition.channel_count,
            tuple(
                (0.0j,) * definition.channel_count
                for _ in range(definition.point_count)
            ),
            "test.matmul",
            "1",
        )


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


def test__init__constructs_default_nominal_evaluators() -> None:
    evaluator = BlockHoppingBlochEvaluator()

    assert isinstance(evaluator.finite_series_evaluator, FiniteFourierSeriesEvaluator)
    assert isinstance(
        evaluator.batched_finite_series_evaluator,
        BatchedFiniteFourierSeriesEvaluator,
    )


def test__init__retains_explicit_batched_evaluator() -> None:
    batched = _BatchedEvaluator()

    evaluator = BlockHoppingBlochEvaluator(
        batched_finite_series_evaluator=batched,
    )

    assert evaluator.batched_finite_series_evaluator is batched


def test__init__uses_scalar_only_when_scalar_owner_is_injected() -> None:
    evaluator = BlockHoppingBlochEvaluator(_Evaluator())

    assert evaluator.batched_finite_series_evaluator is None


def test__init__rejects_nonowned_evaluators() -> None:
    with pytest.raises(TypeError, match="FiniteFourierSeriesEvaluator or None"):
        BlockHoppingBlochEvaluator(cast(Any, object()))
    with pytest.raises(TypeError, match="BatchedFiniteFourierSeriesEvaluator or None"):
        BlockHoppingBlochEvaluator(
            batched_finite_series_evaluator=cast(Any, object()),
        )
