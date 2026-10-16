"""Range-family band evaluator construction tests."""

from typing import Any, cast

import pytest

from projectkoios.physkit.solidstate.tb.bands import (
    BlockHoppingBandEvaluator,
    BlockHoppingRangeBandEvaluator,
)
from projectkoios.physkit.solidstate.tb.truncation import (
    ConjugateClosedBlockHoppingTruncator,
)


def test__init__constructs_or_retains_collaborators() -> None:
    default = BlockHoppingRangeBandEvaluator()
    band_evaluator = BlockHoppingBandEvaluator()
    truncator = ConjugateClosedBlockHoppingTruncator()
    injected = BlockHoppingRangeBandEvaluator(band_evaluator, truncator)

    assert type(default.band_evaluator) is BlockHoppingBandEvaluator
    assert type(default.truncator) is ConjugateClosedBlockHoppingTruncator
    assert injected.band_evaluator is band_evaluator
    assert injected.truncator is truncator


def test__init__rejects_wrong_collaborators() -> None:
    with pytest.raises(TypeError, match="band_evaluator"):
        BlockHoppingRangeBandEvaluator(cast(Any, object()))
    with pytest.raises(TypeError, match="truncator"):
        BlockHoppingRangeBandEvaluator(truncator=cast(Any, object()))
