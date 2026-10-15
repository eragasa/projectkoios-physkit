"""Explicit conjugate-closed block-hopping truncation tests."""

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
    ConjugateClosedBlockHoppingTruncator,
)
from projectkoios.physkit.units import ComplexMatrixQuantity, PhysicalUnit


def _model() -> BlockHoppingModel:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    return BlockHoppingModel(
        "five-term model",
        dimension,
        tuple(
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (representative,)),
                ComplexMatrixQuantity(np.array([[value]]), unit),
            )
            for representative, value in (
                (-2, 5.0),
                (-1, 2.0),
                (0, 3.0),
                (1, 2.0),
                (2, 5.0),
            )
        ),
        unit,
        "common zero",
        "ordered basis",
    )


def test__execute__retains_explicit_terms_metadata_and_omitted_power() -> None:
    source = _model()
    requested = tuple(
        LatticeDisplacement(LatticeDimension.ONE, (value,)) for value in (-1, 0, 1)
    )

    result = ConjugateClosedBlockHoppingTruncator().execute(
        source,
        requested,
        "three-term retained model",
    )

    assert result.source is source
    assert tuple(
        displacement.components for displacement in result.retained_displacements
    ) == ((-1,), (0,), (1,))
    assert all(
        retained is source.terms[index]
        for retained, index in zip(result.truncated.terms, (1, 2, 3), strict=True)
    )
    assert result.truncated.identifier == "three-term retained model"
    assert result.truncated.identifier != source.identifier
    assert result.truncated.dimension is source.dimension
    assert result.truncated.energy_unit == source.energy_unit
    assert result.truncated.energy_reference == source.energy_reference
    assert result.truncated.basis_identifier == source.basis_identifier
    assert tuple(
        displacement.components for displacement in result.omitted_displacements
    ) == ((-2,), (2,))
    assert result.omitted_power.magnitude == 50.0
    assert result.omitted_power.unit == PhysicalUnit("(electron_volt) ** 2")


def test__execute__does_not_require_zero_displacement() -> None:
    source = _model()
    retained = (source.terms[1].displacement, source.terms[3].displacement)

    result = ConjugateClosedBlockHoppingTruncator().execute(
        source,
        retained,
        "paired retained model",
    )

    assert tuple(term.displacement.components for term in result.truncated.terms) == (
        (-1,),
        (1,),
    )


def test__execute__rejects_nonclosed_or_invalid_selection() -> None:
    source = _model()
    action = ConjugateClosedBlockHoppingTruncator()
    negative = source.terms[1].displacement
    zero = source.terms[2].displacement
    positive = source.terms[3].displacement

    with pytest.raises(TypeError, match="source must be BlockHoppingModel"):
        action.execute(cast(Any, object()), (zero,), "retained")
    with pytest.raises(TypeError, match="tuple of LatticeDisplacement"):
        action.execute(source, cast(Any, [zero]), "retained")
    with pytest.raises(ValueError, match="must be nonempty"):
        action.execute(source, (), "retained")
    with pytest.raises(ValueError, match="sorted and unique"):
        action.execute(source, (positive, negative), "retained")
    with pytest.raises(ValueError, match="sorted and unique"):
        action.execute(source, (zero, zero), "retained")
    with pytest.raises(ValueError, match="dimensions must match"):
        action.execute(
            source,
            (LatticeDisplacement(LatticeDimension.TWO, (0, 0)),),
            "retained",
        )
    with pytest.raises(ValueError, match="belong to the source"):
        action.execute(
            source,
            (LatticeDisplacement(LatticeDimension.ONE, (3,)),),
            "retained",
        )
    with pytest.raises(ValueError, match="closed under negation"):
        action.execute(source, (zero, positive), "retained")
    with pytest.raises(TypeError, match="truncated_identifier must be a string"):
        action.execute(source, (zero,), cast(Any, object()))
    with pytest.raises(ValueError, match="must be nonempty"):
        action.execute(source, (zero,), "")
    with pytest.raises(ValueError, match="must differ from source identifier"):
        action.execute(source, (zero,), source.identifier)


def test__execute__does_not_apply_cyclic_nyquist_aliasing() -> None:
    source = _model()
    aliased_terms = source.terms[:-1]
    aliased_source = BlockHoppingModel(
        "centered even-mesh inventory",
        source.dimension,
        aliased_terms,
        source.energy_unit,
        source.energy_reference,
        source.basis_identifier,
    )
    retained = tuple(term.displacement for term in aliased_terms[1:])
    action = ConjugateClosedBlockHoppingTruncator()

    result = action.execute(aliased_source, retained, "non-Nyquist truncation")

    assert tuple(
        displacement.components for displacement in result.omitted_displacements
    ) == ((-2,),)
    with pytest.raises(ValueError, match="closed under negation"):
        action.execute(
            aliased_source,
            (aliased_terms[0].displacement,),
            "invalid Nyquist reinterpretation",
        )
