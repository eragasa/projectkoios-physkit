"""Block-hopping band evaluator construction tests."""

from typing import Any, cast

import pytest

from projectkoios.physkit.solidstate.tb.bands import BlockHoppingBandEvaluator
from projectkoios.physkit.solidstate.tb.bloch import DenseBlockHoppingBlochEvaluator


def test__init__constructs_or_retains_bloch_evaluator() -> None:
    default = BlockHoppingBandEvaluator()
    injected_bloch_evaluator = DenseBlockHoppingBlochEvaluator()
    injected = BlockHoppingBandEvaluator(injected_bloch_evaluator)

    assert type(default.bloch_evaluator) is DenseBlockHoppingBlochEvaluator
    assert injected.bloch_evaluator is injected_bloch_evaluator


def test__init__rejects_wrong_bloch_evaluator() -> None:
    with pytest.raises(TypeError, match="DenseBlockHoppingBlochEvaluator or None"):
        BlockHoppingBandEvaluator(cast(Any, object()))
