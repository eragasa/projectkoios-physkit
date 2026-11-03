"""Invariant tests for correlated matrix Bloch evaluations."""

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
from projectkoios.physkit.solidstate.tb.bloch import DenseBlockHoppingBlochEvaluation
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)


def _model() -> DenseBlockHoppingModel:
    unit = PhysicalUnit("electron_volt")
    return DenseBlockHoppingModel(
        "two-band",
        LatticeDimension.ONE,
        (
            DenseBlockHoppingTerm(
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

    result = DenseBlockHoppingBlochEvaluation(model, coordinates, matrices)

    assert result.model is model
    assert result.fractional_wave_vectors is coordinates
    assert result.bloch_matrices == matrices


def test__init__rejects_wrong_count_shape_unit_and_types() -> None:
    model = _model()
    coordinates = MatrixQuantity(np.array([[0.0], [0.5]]), Unitless())
    valid = ComplexMatrixQuantity(np.eye(2), model.energy_unit)

    with pytest.raises(TypeError, match="model must be DenseBlockHoppingModel"):
        DenseBlockHoppingBlochEvaluation(
            cast(Any, object()), coordinates, (valid, valid)
        )
    with pytest.raises(TypeError, match="tuple of ComplexMatrixQuantity"):
        DenseBlockHoppingBlochEvaluation(model, coordinates, cast(Any, [valid, valid]))
    with pytest.raises(ValueError, match="wave-vector count"):
        DenseBlockHoppingBlochEvaluation(model, coordinates, (valid,))
    with pytest.raises(ValueError, match="basis_size"):
        DenseBlockHoppingBlochEvaluation(
            model,
            coordinates,
            (
                ComplexMatrixQuantity(np.ones((1, 1)), model.energy_unit),
                ComplexMatrixQuantity(np.ones((1, 1)), model.energy_unit),
            ),
        )
    with pytest.raises(ValueError, match="model energy unit"):
        DenseBlockHoppingBlochEvaluation(
            model,
            coordinates,
            (
                ComplexMatrixQuantity(np.eye(2), PhysicalUnit("joule")),
                ComplexMatrixQuantity(np.eye(2), PhysicalUnit("joule")),
            ),
        )
