"""Aggregate block-band diagnostic invariant tests."""

from typing import Any, cast

import pytest

from projectkoios.physkit.solidstate.tb.bands import BlockHoppingBandDiagnostics
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity


def _quantity(value: float, unit: str = "electron_volt") -> ScalarQuantity:
    return ScalarQuantity(value, PhysicalUnit(unit))


def test__init__accepts_componentwise_nonnegative_maxima() -> None:
    diagnostics = BlockHoppingBandDiagnostics(
        maximum_input_hermiticity_absolute_defect=_quantity(0.1),
        maximum_representative_adjustment_absolute_defect=_quantity(0.05),
        maximum_input_matrix_absolute_entry=_quantity(3.0),
        maximum_applied_hermiticity_cutoff=_quantity(0.2),
    )

    assert diagnostics.maximum_input_matrix_absolute_entry.magnitude == 3.0


def test__init__rejects_wrong_quantity_negative_value_and_mixed_units() -> None:
    values = dict(
        maximum_input_hermiticity_absolute_defect=_quantity(0.1),
        maximum_representative_adjustment_absolute_defect=_quantity(0.05),
        maximum_input_matrix_absolute_entry=_quantity(3.0),
        maximum_applied_hermiticity_cutoff=_quantity(0.2),
    )
    with pytest.raises(TypeError, match="ScalarQuantity"):
        BlockHoppingBandDiagnostics(
            **{
                **values,
                "maximum_input_matrix_absolute_entry": cast(Any, 3.0),
            }
        )
    with pytest.raises(ValueError, match="nonnegative"):
        BlockHoppingBandDiagnostics(
            **{
                **values,
                "maximum_representative_adjustment_absolute_defect": _quantity(-0.1),
            }
        )
    with pytest.raises(ValueError, match="common unit"):
        BlockHoppingBandDiagnostics(
            **{
                **values,
                "maximum_applied_hermiticity_cutoff": _quantity(0.2, "joule"),
            }
        )
    with pytest.raises(ValueError, match="must not exceed"):
        BlockHoppingBandDiagnostics(
            **{
                **values,
                "maximum_applied_hermiticity_cutoff": _quantity(0.05),
            }
        )
