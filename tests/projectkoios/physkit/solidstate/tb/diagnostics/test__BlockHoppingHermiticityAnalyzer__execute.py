"""Raw Hermiticity diagnostics for matrix-valued hopping coefficients."""

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
from projectkoios.physkit.solidstate.tb.diagnostics import (
    BlockHoppingHermiticityAnalyzer,
)
from projectkoios.physkit.units import ComplexMatrixQuantity, PhysicalUnit


def _model(
    entries: tuple[tuple[tuple[int, int], np.ndarray], ...],
) -> BlockHoppingModel:
    dimension = LatticeDimension.TWO
    unit = PhysicalUnit("electron_volt")
    return BlockHoppingModel(
        "block diagnostic",
        dimension,
        tuple(
            BlockHoppingTerm(
                LatticeDisplacement(dimension, displacement),
                ComplexMatrixQuantity(matrix, unit),
            )
            for displacement, matrix in entries
        ),
        unit,
        "zero",
        "ordered basis",
    )


def test__execute__compares_pairs_once_and_reports_missing_targets() -> None:
    negative = np.array([[1.0, 2.0j], [3.0, 4.0]], dtype=np.complex128)
    positive = negative.conj().T.copy()
    positive[0, 0] += 0.5
    model = _model(
        (
            ((-1, 0), negative),
            ((0, 0), np.diag([2.0, 3.0]).astype(np.complex128)),
            ((0, 1), np.eye(2, dtype=np.complex128)),
            ((1, 0), positive),
        )
    )

    result = BlockHoppingHermiticityAnalyzer().execute(model)

    assert result.model is model
    assert tuple(
        (
            comparison.source.components,
            comparison.target.components,
            None
            if comparison.frobenius_defect is None
            else comparison.frobenius_defect.magnitude,
        )
        for comparison in result.comparisons
    ) == (
        ((-1, 0), (1, 0), 0.5),
        ((0, 0), (0, 0), 0.0),
        ((0, 1), (0, -1), None),
    )
    assert tuple(target.components for target in result.missing_targets) == ((0, -1),)
    assert result.maximum_frobenius_defect is not None
    assert result.maximum_frobenius_defect.magnitude == 0.5
    assert result.maximum_frobenius_defect.unit == model.energy_unit


def test__execute__reports_no_maximum_when_every_target_is_missing() -> None:
    result = BlockHoppingHermiticityAnalyzer().execute(
        _model((((1, 0), np.eye(2, dtype=np.complex128)),))
    )

    assert result.maximum_frobenius_defect is None
    assert tuple(target.components for target in result.missing_targets) == ((-1, 0),)


def test__execute__does_not_apply_cyclic_nyquist_aliasing() -> None:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    model = BlockHoppingModel(
        "centered even-mesh inventory",
        dimension,
        (
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (-2,)),
                ComplexMatrixQuantity(np.eye(2, dtype=np.complex128), unit),
            ),
        ),
        unit,
        "zero",
        "basis",
    )

    result = BlockHoppingHermiticityAnalyzer().execute(model)

    assert tuple(target.components for target in result.missing_targets) == ((2,),)
    assert result.maximum_frobenius_defect is None


def test__execute__rejects_wrong_model_and_nonfinite_defect_arithmetic() -> None:
    analyzer = BlockHoppingHermiticityAnalyzer()
    with pytest.raises(TypeError, match="model must be BlockHoppingModel"):
        analyzer.execute(cast(Any, object()))

    negative = np.full((1, 1), -1.0e308, dtype=np.complex128)
    positive = np.full((1, 1), 1.0e308, dtype=np.complex128)
    with pytest.raises(FloatingPointError, match="Frobenius"):
        analyzer.execute(_model((((-1, 0), negative), ((1, 0), positive))))
