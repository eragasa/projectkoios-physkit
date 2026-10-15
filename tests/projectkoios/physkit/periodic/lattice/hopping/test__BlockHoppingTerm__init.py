"""Construction tests for matrix-valued hopping terms."""

from dataclasses import FrozenInstanceError
from typing import Any, cast

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import BlockHoppingTerm
from projectkoios.physkit.units import ComplexMatrixQuantity, PhysicalUnit


def test__init__retains_displacement_and_immutable_coefficient() -> None:
    displacement = LatticeDisplacement(LatticeDimension.TWO, (-1, 2))
    coefficient = ComplexMatrixQuantity(
        np.array([[1.0, 2.0j], [-2.0j, 3.0]], dtype=np.complex128),
        PhysicalUnit("electron_volt"),
    )

    term = BlockHoppingTerm(displacement, coefficient)

    assert term.displacement is displacement
    assert term.coefficient is coefficient
    with pytest.raises(FrozenInstanceError):
        term.displacement = LatticeDisplacement(  # type: ignore[misc]
            LatticeDimension.TWO,
            (0, 0),
        )
    with pytest.raises(ValueError, match="read-only"):
        term.coefficient.magnitude[0, 0] = 0.0


def test__init__rejects_nonowned_inputs() -> None:
    coefficient = ComplexMatrixQuantity(
        np.eye(2, dtype=np.complex128),
        PhysicalUnit("electron_volt"),
    )

    with pytest.raises(TypeError, match="displacement must be LatticeDisplacement"):
        BlockHoppingTerm(cast(Any, object()), coefficient)
    with pytest.raises(TypeError, match="coefficient must be ComplexMatrixQuantity"):
        BlockHoppingTerm(
            LatticeDisplacement(LatticeDimension.ONE, (0,)),
            cast(Any, object()),
        )


@pytest.mark.parametrize(
    "values",
    (
        np.empty((0, 0), dtype=np.complex128),
        np.zeros((2, 3), dtype=np.complex128),
    ),
)
def test__init__rejects_empty_or_nonsquare_coefficients(values: np.ndarray) -> None:
    coefficient = ComplexMatrixQuantity(values, PhysicalUnit("electron_volt"))

    with pytest.raises(ValueError, match="nonempty square"):
        BlockHoppingTerm(
            LatticeDisplacement(LatticeDimension.ONE, (0,)),
            coefficient,
        )
