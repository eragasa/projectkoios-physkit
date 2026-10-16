"""Explicit range-family selection tests."""

from typing import Any, cast

import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.solidstate.tb.bands import BlockHoppingRangeBandSelection


def _displacements(*values: int) -> tuple[LatticeDisplacement, ...]:
    return tuple(
        LatticeDisplacement(LatticeDimension.ONE, (value,)) for value in values
    )


def test__init__retains_explicit_inventory_and_identifier() -> None:
    retained = _displacements(-1, 0, 1)

    selection = BlockHoppingRangeBandSelection(retained, "range one")

    assert selection.retained_displacements is retained
    assert selection.truncated_identifier == "range one"


def test__init__rejects_invalid_inventory_or_identifier() -> None:
    with pytest.raises(TypeError, match="tuple of LatticeDisplacement"):
        BlockHoppingRangeBandSelection(cast(Any, [_displacements(0)[0]]), "range")
    with pytest.raises(ValueError, match="nonempty"):
        BlockHoppingRangeBandSelection((), "range")
    with pytest.raises(ValueError, match="sorted and unique"):
        BlockHoppingRangeBandSelection(_displacements(0, -1, 1), "range")
    with pytest.raises(ValueError, match="sorted and unique"):
        BlockHoppingRangeBandSelection(_displacements(0, 0), "range")
    with pytest.raises(ValueError, match="closed under negation"):
        BlockHoppingRangeBandSelection(_displacements(0, 1), "range")
    with pytest.raises(ValueError, match="dimensions must match"):
        BlockHoppingRangeBandSelection(
            (
                LatticeDisplacement(LatticeDimension.ONE, (0,)),
                LatticeDisplacement(LatticeDimension.TWO, (0, 0)),
            ),
            "range",
        )
    with pytest.raises(TypeError, match="must be a string"):
        BlockHoppingRangeBandSelection(_displacements(0), cast(Any, object()))
    with pytest.raises(ValueError, match="must be nonempty"):
        BlockHoppingRangeBandSelection(_displacements(0), "")
