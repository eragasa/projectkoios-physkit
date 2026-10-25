"""Raw unit-aware block Frobenius-power analysis tests."""

from typing import Any, cast

import numpy as np
import pytest

from projectkoios.numerics.linalg.numpy.frobenius import NumpyFrobeniusPower
from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    BlockHoppingModel,
    BlockHoppingTerm,
)
from projectkoios.physkit.solidstate.tb import diagnostics as diagnostics_module
from projectkoios.physkit.solidstate.tb.diagnostics import BlockHoppingPowerAnalyzer
from projectkoios.physkit.units import ComplexMatrixQuantity, PhysicalUnit, Unitless


def _model(unit: PhysicalUnit | Unitless) -> BlockHoppingModel:
    dimension = LatticeDimension.ONE
    return BlockHoppingModel(
        "power model",
        dimension,
        (
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (0,)),
                ComplexMatrixQuantity(
                    np.array([[3.0, 4.0], [0.0, 0.0]]),
                    unit,
                ),
            ),
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (1,)),
                ComplexMatrixQuantity(
                    np.array([[1.0j, 2.0], [2.0, 0.0]]),
                    unit,
                ),
            ),
        ),
        unit,
        "zero",
        "basis",
    )


def test__execute__returns_raw_powers_in_squared_energy_units() -> None:
    model = _model(PhysicalUnit("electron_volt"))

    result = BlockHoppingPowerAnalyzer().execute(model)

    assert result.model is model
    assert tuple(power.frobenius_power.magnitude for power in result.powers) == (
        25.0,
        9.0,
    )
    assert all(
        power.frobenius_power.unit == PhysicalUnit("(electron_volt) ** 2")
        for power in result.powers
    )
    assert result.total_power.magnitude == 34.0
    assert result.total_power.unit == PhysicalUnit("(electron_volt) ** 2")


def test__execute__delegates_exact_dense_matrices_to_numerics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[np.ndarray] = []

    def recording_frobenius_power(matrix: np.ndarray, /) -> NumpyFrobeniusPower:
        calls.append(matrix)
        return NumpyFrobeniusPower(matrix)

    monkeypatch.setattr(
        diagnostics_module,
        "NumpyFrobeniusPower",
        recording_frobenius_power,
    )

    model = _model(PhysicalUnit("electron_volt"))
    BlockHoppingPowerAnalyzer().execute(model)

    expected = tuple(term.coefficient.magnitude for term in model.terms) * 2
    assert len(calls) == len(expected)
    assert all(
        actual is retained for actual, retained in zip(calls, expected, strict=True)
    )


def test__execute__preserves_explicit_unitless_power() -> None:
    result = BlockHoppingPowerAnalyzer().execute(_model(Unitless()))

    assert all(
        isinstance(power.frobenius_power.unit, Unitless) for power in result.powers
    )
    assert isinstance(result.total_power.unit, Unitless)


def test__execute__rejects_wrong_model_and_nonfinite_power() -> None:
    analyzer = BlockHoppingPowerAnalyzer()
    with pytest.raises(TypeError, match="model must be BlockHoppingModel"):
        analyzer.execute(cast(Any, object()))

    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    model = BlockHoppingModel(
        "overflow",
        dimension,
        (
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (0,)),
                ComplexMatrixQuantity(np.array([[1.0e308]]), unit),
            ),
        ),
        unit,
        "zero",
        "basis",
    )
    with pytest.raises(FloatingPointError, match="Frobenius power"):
        analyzer.execute(model)
