"""Block-hopping band result invariant tests."""

from dataclasses import replace
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
from projectkoios.physkit.solidstate.tb.bands import (
    BlockHoppingBandEvaluation,
    BlockHoppingBandEvaluator,
    BlockHoppingBandTolerance,
)
from projectkoios.physkit.solidstate.tb.bloch import BlockHoppingBlochEvaluator
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
)


def _canonical() -> BlockHoppingBandEvaluation:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    model = BlockHoppingModel(
        "Hermitian two-band model",
        dimension,
        (
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (0,)),
                ComplexMatrixQuantity(np.diag([1.0, 2.0]), unit),
            ),
        ),
        unit,
        "zero",
        "basis",
    )
    return BlockHoppingBandEvaluator().evaluate(
        model,
        MatrixQuantity(np.array([[0.0], [0.5]]), Unitless()),
        BlockHoppingBandTolerance(ScalarQuantity(0.0, unit), 0.0),
    )


def test__init__accepts_canonical_correlated_result() -> None:
    canonical = _canonical()

    rebuilt = BlockHoppingBandEvaluation(
        canonical.bloch_evaluation,
        canonical.band_energies,
        canonical.tolerance,
        canonical.point_diagnostics,
        canonical.diagnostics,
        canonical.backend_identifier,
        canonical.backend_version,
    )

    assert rebuilt.bloch_evaluation is canonical.bloch_evaluation
    assert rebuilt.model is canonical.model
    assert rebuilt.fractional_wave_vectors is canonical.fractional_wave_vectors
    assert rebuilt.band_energies is canonical.band_energies


def test__init__rejects_wrong_bloch_evaluation_and_band_contracts() -> None:
    canonical = _canonical()
    unit = canonical.model.energy_unit

    with pytest.raises(TypeError, match="bloch_evaluation"):
        replace(canonical, bloch_evaluation=cast(Any, object()))
    shorter_bloch_evaluation = BlockHoppingBlochEvaluator().evaluate(
        canonical.model,
        MatrixQuantity(np.array([[0.0]]), Unitless()),
    )
    with pytest.raises(ValueError, match="shape must be"):
        replace(canonical, bloch_evaluation=shorter_bloch_evaluation)
    with pytest.raises(ValueError, match="shape must be"):
        replace(
            canonical,
            band_energies=MatrixQuantity(np.array([[1.0], [1.0]]), unit),
        )
    with pytest.raises(ValueError, match="model energy unit"):
        replace(
            canonical,
            band_energies=MatrixQuantity(
                canonical.band_energies.magnitude,
                PhysicalUnit("joule"),
            ),
        )
    with pytest.raises(ValueError, match="rows must be ascending"):
        replace(
            canonical,
            band_energies=MatrixQuantity(
                np.array([[2.0, 1.0], [1.0, 2.0]]),
                unit,
            ),
        )


def test__init__rejects_uncorrelated_tolerance_and_diagnostics() -> None:
    canonical = _canonical()
    unit = canonical.model.energy_unit

    with pytest.raises(ValueError, match="absolute tolerance"):
        replace(
            canonical,
            tolerance=BlockHoppingBandTolerance(
                ScalarQuantity(0.0, PhysicalUnit("joule")),
                0.0,
            ),
        )
    with pytest.raises(ValueError, match="Bloch-matrix count"):
        replace(canonical, point_diagnostics=())
    with pytest.raises(ValueError, match="canonical point indices"):
        replace(
            canonical,
            point_diagnostics=(
                replace(canonical.point_diagnostics[0], point_index=1),
                canonical.point_diagnostics[1],
            ),
        )
    with pytest.raises(ValueError, match="componentwise point maxima"):
        replace(canonical, diagnostics=None)
    assert canonical.diagnostics is not None
    wrong_aggregate = replace(
        canonical.diagnostics,
        maximum_input_matrix_absolute_entry=ScalarQuantity(3.0, unit),
    )
    with pytest.raises(ValueError, match="componentwise point maxima"):
        replace(canonical, diagnostics=wrong_aggregate)
    wrong_cutoff_points = (
        replace(
            canonical.point_diagnostics[0],
            applied_hermiticity_cutoff=ScalarQuantity(0.1, unit),
        ),
        canonical.point_diagnostics[1],
    )
    wrong_cutoff_aggregate = replace(
        canonical.diagnostics,
        maximum_applied_hermiticity_cutoff=ScalarQuantity(0.1, unit),
    )
    with pytest.raises(ValueError, match="applied cutoffs must match"):
        replace(
            canonical,
            point_diagnostics=wrong_cutoff_points,
            diagnostics=wrong_cutoff_aggregate,
        )


def test__init__rejects_missing_backend_evidence() -> None:
    canonical = _canonical()

    with pytest.raises(TypeError, match="backend_identifier must be a string"):
        replace(canonical, backend_identifier=cast(Any, 1))
    with pytest.raises(ValueError, match="backend_identifier must be nonempty"):
        replace(canonical, backend_identifier="")
    with pytest.raises(TypeError, match="backend_version must be a string"):
        replace(canonical, backend_version=cast(Any, 1))
    with pytest.raises(ValueError, match="backend_version must be nonempty"):
        replace(canonical, backend_version="")
