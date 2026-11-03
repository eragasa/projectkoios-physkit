"""Correlation tests for block-hopping Hermiticity results."""

from typing import Any, cast

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    DenseBlockHoppingModel,
    DenseBlockHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.diagnostics import (
    BlockHoppingCoefficientComparison,
    BlockHoppingHermiticityAnalyzer,
    BlockHoppingHermiticityResult,
)
from projectkoios.physkit.units import ComplexMatrixQuantity, PhysicalUnit


def _model() -> DenseBlockHoppingModel:
    unit = PhysicalUnit("electron_volt")
    dimension = LatticeDimension.ONE
    term = DenseBlockHoppingTerm(
        LatticeDisplacement(dimension, (0,)),
        ComplexMatrixQuantity(np.eye(2), unit),
    )
    return DenseBlockHoppingModel("model", dimension, (term,), unit, "zero", "basis")


def test__init__accepts_only_canonical_comparisons() -> None:
    canonical = BlockHoppingHermiticityAnalyzer().execute(_model())

    rebuilt = BlockHoppingHermiticityResult(
        canonical.model,
        canonical.comparisons,
    )

    assert rebuilt.comparisons == canonical.comparisons


def test__init__rejects_wrong_types_and_noncanonical_values() -> None:
    model = _model()
    canonical = BlockHoppingHermiticityAnalyzer().execute(model)

    with pytest.raises(TypeError, match="model must be DenseBlockHoppingModel"):
        BlockHoppingHermiticityResult(cast(Any, object()), canonical.comparisons)
    with pytest.raises(TypeError, match="BlockHoppingCoefficientComparison"):
        BlockHoppingHermiticityResult(model, cast(Any, (object(),)))
    fabricated = BlockHoppingCoefficientComparison(
        canonical.comparisons[0].source,
        canonical.comparisons[0].target,
        None,
    )
    with pytest.raises(ValueError, match="canonical block comparisons"):
        BlockHoppingHermiticityResult(model, (fabricated,))
