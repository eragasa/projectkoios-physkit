"""Invariant tests for correlated matrix Bloch evaluations."""

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
from projectkoios.physkit.solidstate.tb.bloch import BlockHoppingBlochEvaluation
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)


def _model() -> BlockHoppingModel:
    unit = PhysicalUnit("electron_volt")
    return BlockHoppingModel(
        "two-band",
        LatticeDimension.ONE,
        (
            BlockHoppingTerm(
                LatticeDisplacement(LatticeDimension.ONE, (0,)),
                ComplexMatrixQuantity(np.eye(2), unit),
            ),
        ),
        unit,
        "zero",
        "basis",
    )


def test__init__retains_correlated_model_coordinates_and_matrices() -> None:
    model = _model()
    coordinates = MatrixQuantity(np.array([[0.0], [0.5]]), Unitless())
    matrices = tuple(
        ComplexMatrixQuantity(np.eye(2) * value, model.energy_unit)
        for value in (1.0, 2.0)
    )

    result = BlockHoppingBlochEvaluation(model, coordinates, matrices)

    assert result.model is model
    assert result.fractional_wave_vectors is coordinates
    assert result.bloch_matrices == matrices


def test__init__rejects_wrong_count_shape_unit_and_types() -> None:
    model = _model()
    coordinates = MatrixQuantity(np.array([[0.0], [0.5]]), Unitless())
    valid = ComplexMatrixQuantity(np.eye(2), model.energy_unit)

    with pytest.raises(TypeError, match="model must be BlockHoppingModel"):
        BlockHoppingBlochEvaluation(cast(Any, object()), coordinates, (valid, valid))
    with pytest.raises(TypeError, match="tuple of ComplexMatrixQuantity"):
        BlockHoppingBlochEvaluation(model, coordinates, cast(Any, [valid, valid]))
    with pytest.raises(ValueError, match="wave-vector count"):
        BlockHoppingBlochEvaluation(model, coordinates, (valid,))
    with pytest.raises(ValueError, match="basis_size"):
        BlockHoppingBlochEvaluation(
            model,
            coordinates,
            (
                ComplexMatrixQuantity(np.ones((1, 1)), model.energy_unit),
                ComplexMatrixQuantity(np.ones((1, 1)), model.energy_unit),
            ),
        )
    with pytest.raises(ValueError, match="model energy unit"):
        BlockHoppingBlochEvaluation(
            model,
            coordinates,
            (
                ComplexMatrixQuantity(np.eye(2), PhysicalUnit("joule")),
                ComplexMatrixQuantity(np.eye(2), PhysicalUnit("joule")),
            ),
        )
