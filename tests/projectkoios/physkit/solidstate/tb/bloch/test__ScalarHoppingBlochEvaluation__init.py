"""Construction tests for scalar-hopping Bloch evaluations."""

from __future__ import annotations

from typing import Any, cast

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    ScalarHoppingModel,
    ScalarHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.bloch import ScalarHoppingBlochEvaluation
from projectkoios.physkit.units import (
    ComplexVectorQuantity,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)

pytestmark = pytest.mark.unit


def _model() -> ScalarHoppingModel:
    return ScalarHoppingModel(
        identifier="two-dimensional scalar model",
        dimension=LatticeDimension.TWO,
        terms=(
            ScalarHoppingTerm(
                LatticeDisplacement(LatticeDimension.TWO, (0, 0)),
                1.0,
                0.0,
            ),
        ),
        energy_unit=PhysicalUnit("electron_volt"),
        energy_reference="onsite zero",
        basis_identifier="one scalar orbital",
    )


def test__init__retains_model_coordinates_and_energy_unit() -> None:
    model = _model()
    fractional_wave_vectors = MatrixQuantity(
        np.array([[0.0, 0.0], [0.25, -0.5]]),
        Unitless(),
    )
    bloch_values = ComplexVectorQuantity(
        np.array([1.0 + 0.0j, 2.0 - 0.5j]),
        model.energy_unit,
    )

    evaluation = ScalarHoppingBlochEvaluation(
        model=model,
        fractional_wave_vectors=fractional_wave_vectors,
        bloch_values=bloch_values,
    )

    assert evaluation.model is model
    assert evaluation.fractional_wave_vectors is fractional_wave_vectors
    assert evaluation.bloch_values is bloch_values


def test__init__requires_fractional_wave_vector_quantity() -> None:
    model = _model()

    with pytest.raises(
        TypeError, match="fractional_wave_vectors must be MatrixQuantity"
    ):
        ScalarHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=cast(Any, object()),
            bloch_values=ComplexVectorQuantity(
                np.zeros(1, dtype=np.complex128),
                model.energy_unit,
            ),
        )


def test__init__requires_unitless_fractional_wave_vectors() -> None:
    model = _model()

    with pytest.raises(ValueError, match="fractional_wave_vectors must be unitless"):
        ScalarHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=MatrixQuantity(
                np.zeros((1, 2)),
                PhysicalUnit("1 / meter"),
            ),
            bloch_values=ComplexVectorQuantity(
                np.zeros(1, dtype=np.complex128),
                model.energy_unit,
            ),
        )


def test__init__requires_wave_vector_width_to_match_model_dimension() -> None:
    model = _model()

    with pytest.raises(ValueError, match=r"shape \(point_count, 2\)"):
        ScalarHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=MatrixQuantity(
                np.zeros((1, 3)),
                Unitless(),
            ),
            bloch_values=ComplexVectorQuantity(
                np.zeros(1, dtype=np.complex128),
                model.energy_unit,
            ),
        )


def test__init__requires_complex_vector_bloch_values() -> None:
    model = _model()

    with pytest.raises(TypeError, match="bloch_values must be ComplexVectorQuantity"):
        ScalarHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=MatrixQuantity(
                np.zeros((1, 2)),
                Unitless(),
            ),
            bloch_values=cast(Any, object()),
        )


def test__init__requires_bloch_values_to_use_model_energy_unit() -> None:
    model = _model()

    with pytest.raises(ValueError, match="model energy unit"):
        ScalarHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=MatrixQuantity(
                np.zeros((1, 2)),
                Unitless(),
            ),
            bloch_values=ComplexVectorQuantity(
                np.zeros(1, dtype=np.complex128),
                PhysicalUnit("joule"),
            ),
        )


def test__init__requires_one_bloch_value_per_wave_vector() -> None:
    model = _model()

    with pytest.raises(ValueError, match="wave-vector count"):
        ScalarHoppingBlochEvaluation(
            model=model,
            fractional_wave_vectors=MatrixQuantity(
                np.zeros((2, 2)),
                Unitless(),
            ),
            bloch_values=ComplexVectorQuantity(
                np.zeros(1, dtype=np.complex128),
                model.energy_unit,
            ),
        )


def test__init__rejects_nonowned_model_type() -> None:
    model = _model()

    with pytest.raises(TypeError, match="model must be ScalarHoppingModel"):
        ScalarHoppingBlochEvaluation(
            model=cast(Any, object()),
            fractional_wave_vectors=MatrixQuantity(
                np.zeros((1, 2)),
                Unitless(),
            ),
            bloch_values=ComplexVectorQuantity(
                np.zeros(1, dtype=np.complex128),
                model.energy_unit,
            ),
        )
