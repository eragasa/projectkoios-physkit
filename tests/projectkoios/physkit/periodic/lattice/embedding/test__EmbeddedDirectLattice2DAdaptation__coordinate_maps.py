"""Failure and storage checks for embedded-lattice coordinate maps."""

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.embedding import (
    EmbeddedDirectLattice2D,
    EmbeddedDirectLattice2DAdaptationRequest,
    EmbeddedDirectLattice2DAdapter,
)
from projectkoios.physkit.units import MatrixQuantity, Unitless


def _adaptation():
    return EmbeddedDirectLattice2DAdapter().action(
        request=EmbeddedDirectLattice2DAdaptationRequest(
            embedded_lattice=EmbeddedDirectLattice2D(
                ambient_basis=MatrixQuantity(
                    np.array(
                        [
                            [2.0, 0.4],
                            [0.3, 1.7],
                            [1.1, -0.2],
                        ]
                    ),
                    Unitless(),
                )
            )
        )
    )


def test__coordinate_maps__reject_wrong_shape_type_and_nonfinite_values() -> None:
    result = _adaptation()

    with pytest.raises(ValueError, match="final dimension 2"):
        result.fractional_to_ambient(np.ones((2, 3)))
    with pytest.raises(ValueError, match="final dimension 3"):
        result.project_ambient_to_fractional(np.ones((2, 2)))
    with pytest.raises(TypeError, match="numpy.ndarray"):
        result.intrinsic_to_ambient([0.0, 0.0])  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="finite values"):
        result.project_ambient_to_intrinsic(np.array([0.0, np.nan, 0.0]))


def test__coordinate_maps__retain_immutable_correlated_geometry() -> None:
    result = _adaptation()

    for quantity in (
        result.intrinsic_basis,
        result.induced_metric,
        result.inverse_induced_metric,
        result.dual_reciprocal_basis,
        result.orthonormal_frame,
        result.ambient_projector,
    ):
        assert quantity.magnitude.flags.writeable is False

    assert result.intrinsic_metric.primitive_basis.flags.writeable is False
