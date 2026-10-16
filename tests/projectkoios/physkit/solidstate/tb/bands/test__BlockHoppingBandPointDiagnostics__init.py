"""Per-point block-band diagnostic invariant tests."""

from typing import Any, cast

import pytest

from projectkoios.physkit.solidstate.tb.bands import (
    BlockHoppingBandPointDiagnostics,
)
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity


def _quantity(value: float, unit: str = "electron_volt") -> ScalarQuantity:
    return ScalarQuantity(value, PhysicalUnit(unit))


def test__init__accepts_nonnegative_common_unit_diagnostics() -> None:
    diagnostic = BlockHoppingBandPointDiagnostics(
        point_index=2,
        input_hermiticity_maximum_absolute_defect=_quantity(0.1),
        representative_adjustment_maximum_absolute_defect=_quantity(0.05),
        input_matrix_maximum_absolute_entry=_quantity(3.0),
        applied_hermiticity_cutoff=_quantity(0.2),
    )

    assert diagnostic.point_index == 2
    assert diagnostic.input_hermiticity_maximum_absolute_defect.magnitude == 0.1


def test__init__rejects_invalid_index_quantity_unit_and_cutoff() -> None:
    values = dict(
        input_hermiticity_maximum_absolute_defect=_quantity(0.1),
        representative_adjustment_maximum_absolute_defect=_quantity(0.05),
        input_matrix_maximum_absolute_entry=_quantity(3.0),
        applied_hermiticity_cutoff=_quantity(0.2),
    )
    with pytest.raises(TypeError, match="point_index"):
        BlockHoppingBandPointDiagnostics(point_index=cast(Any, 0.0), **values)
    with pytest.raises(ValueError, match="point_index"):
        BlockHoppingBandPointDiagnostics(point_index=-1, **values)
    with pytest.raises(TypeError, match="ScalarQuantity"):
        BlockHoppingBandPointDiagnostics(
            point_index=0,
            **{
                **values,
                "input_matrix_maximum_absolute_entry": cast(Any, 3.0),
            },
        )
    with pytest.raises(ValueError, match="nonnegative"):
        BlockHoppingBandPointDiagnostics(
            point_index=0,
            **{
                **values,
                "representative_adjustment_maximum_absolute_defect": _quantity(-0.1),
            },
        )
    with pytest.raises(ValueError, match="common unit"):
        BlockHoppingBandPointDiagnostics(
            point_index=0,
            **{
                **values,
                "applied_hermiticity_cutoff": _quantity(0.2, "joule"),
            },
        )
    with pytest.raises(ValueError, match="must not exceed"):
        BlockHoppingBandPointDiagnostics(
            point_index=0,
            **{
                **values,
                "applied_hermiticity_cutoff": _quantity(0.05),
            },
        )
