"""Invariant tests for explicit block-hopping truncation results."""

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
from projectkoios.physkit.solidstate.tb.truncation import (
    BlockHoppingTruncationResult,
    ConjugateClosedBlockHoppingTruncator,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
)


def _source() -> BlockHoppingModel:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    return BlockHoppingModel(
        "source",
        dimension,
        tuple(
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (representative,)),
                ComplexMatrixQuantity(np.array([[float(representative + 3)]]), unit),
            )
            for representative in (-1, 0, 1)
        ),
        unit,
        "zero",
        "basis",
    )


def test__init__accepts_canonical_correlated_result() -> None:
    source = _source()
    canonical = ConjugateClosedBlockHoppingTruncator().execute(
        source,
        (source.terms[1].displacement,),
        "retained",
    )

    rebuilt = BlockHoppingTruncationResult(
        canonical.source,
        canonical.retained_displacements,
        canonical.truncated,
        canonical.omitted_power,
    )

    assert rebuilt.truncated is canonical.truncated
    assert rebuilt.omitted_power == canonical.omitted_power


def test__init__rejects_wrong_truncated_model_metadata_terms_and_power() -> None:
    source = _source()
    canonical = ConjugateClosedBlockHoppingTruncator().execute(
        source,
        (source.terms[1].displacement,),
        "retained",
    )
    unchanged_identity = BlockHoppingModel(
        source.identifier,
        source.dimension,
        canonical.truncated.terms,
        source.energy_unit,
        source.energy_reference,
        source.basis_identifier,
    )
    wrong_metadata = BlockHoppingModel(
        "wrong metadata",
        source.dimension,
        canonical.truncated.terms,
        source.energy_unit,
        "different zero",
        source.basis_identifier,
    )
    equivalent_new_term = BlockHoppingTerm(
        source.terms[1].displacement,
        ComplexMatrixQuantity(
            source.terms[1].coefficient.magnitude,
            source.energy_unit,
        ),
    )
    wrong_term_identity = BlockHoppingModel(
        "wrong term identity",
        source.dimension,
        (equivalent_new_term,),
        source.energy_unit,
        source.energy_reference,
        source.basis_identifier,
    )

    with pytest.raises(TypeError, match="truncated must be BlockHoppingModel"):
        BlockHoppingTruncationResult(
            source,
            canonical.retained_displacements,
            cast(Any, object()),
            canonical.omitted_power,
        )
    with pytest.raises(ValueError, match="identifier must differ"):
        BlockHoppingTruncationResult(
            source,
            canonical.retained_displacements,
            unchanged_identity,
            canonical.omitted_power,
        )
    with pytest.raises(ValueError, match="metadata must match"):
        BlockHoppingTruncationResult(
            source,
            canonical.retained_displacements,
            wrong_metadata,
            canonical.omitted_power,
        )
    with pytest.raises(ValueError, match="selected source terms"):
        BlockHoppingTruncationResult(
            source,
            canonical.retained_displacements,
            wrong_term_identity,
            canonical.omitted_power,
        )
    with pytest.raises(TypeError, match="omitted_power must be ScalarQuantity"):
        BlockHoppingTruncationResult(
            source,
            canonical.retained_displacements,
            canonical.truncated,
            cast(Any, 13.0),
        )
    with pytest.raises(ValueError, match="unselected source terms"):
        BlockHoppingTruncationResult(
            source,
            canonical.retained_displacements,
            canonical.truncated,
            ScalarQuantity(13.0, PhysicalUnit("(electron_volt) ** 2")),
        )
