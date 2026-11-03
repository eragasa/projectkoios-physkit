"""Construction tests for sparse-block hopping models."""

from typing import Any, cast

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    SparseBlockHoppingModel,
    SparseBlockHoppingTerm,
)
from projectkoios.physkit.units import (
    ComplexSparseMatrixQuantity,
    PhysicalUnit,
)


def _term(
    displacement: int,
    diagonal: tuple[complex, ...],
    unit: PhysicalUnit,
) -> SparseBlockHoppingTerm:
    return SparseBlockHoppingTerm(
        LatticeDisplacement(LatticeDimension.ONE, (displacement,)),
        ComplexSparseMatrixQuantity.from_csr(
            sparse.diags(np.asarray(diagonal, dtype=np.complex128), format="csr"),
            unit,
        ),
    )


def test__init__ordered_sparse_blocks__retains_identity_and_shape() -> None:
    """A sparse model retains canonical ordered blocks and scientific identities."""
    unit = PhysicalUnit("electron_volt")
    terms = (
        _term(-1, (1.0 - 1.0j, 0.0j), unit),
        _term(0, (2.0 + 0.0j, 3.0 + 0.0j), unit),
        _term(1, (1.0 + 1.0j, 0.0j), unit),
    )

    model = SparseBlockHoppingModel(
        identifier="sparse-two-state-chain",
        dimension=LatticeDimension.ONE,
        terms=terms,
        energy_unit=unit,
        energy_reference="common-zero",
        basis_identifier="ordered basis: a, b",
    )

    assert model.terms == terms
    assert model.basis_size == 2
    assert model.storage_kind == "sparse"
    assert model.energy_unit is unit


def test__init__mixed_shape_unit_dimension_or_order__rejects() -> None:
    """Sparse model correlation fails closed for incompatible term inventories."""
    unit = PhysicalUnit("electron_volt")
    valid = _term(0, (1.0 + 0.0j, 2.0 + 0.0j), unit)
    wrong_shape = _term(1, (1.0 + 0.0j,), unit)
    wrong_unit = _term(1, (1.0 + 0.0j, 2.0 + 0.0j), PhysicalUnit("joule"))
    wrong_dimension = SparseBlockHoppingTerm(
        LatticeDisplacement(LatticeDimension.TWO, (1, 0)),
        valid.coefficient,
    )

    for term, message in (
        (wrong_shape, "common square shape"),
        (wrong_unit, "units must match"),
        (wrong_dimension, "dimensions must match"),
    ):
        with pytest.raises(ValueError, match=message):
            SparseBlockHoppingModel(
                "invalid",
                LatticeDimension.ONE,
                (valid, term),
                unit,
                "zero",
                "basis",
            )
    with pytest.raises(ValueError, match="sorted and unique"):
        SparseBlockHoppingModel(
            "invalid",
            LatticeDimension.ONE,
            (_term(1, (1.0 + 0.0j, 2.0 + 0.0j), unit), valid),
            unit,
            "zero",
            "basis",
        )
    with pytest.raises(TypeError, match="SparseBlockHoppingTerm"):
        SparseBlockHoppingModel(
            "invalid",
            LatticeDimension.ONE,
            (cast(Any, object()),),
            unit,
            "zero",
            "basis",
        )
