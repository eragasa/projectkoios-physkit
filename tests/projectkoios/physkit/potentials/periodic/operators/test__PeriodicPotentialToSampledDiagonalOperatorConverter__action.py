"""Sampled periodic-potential operator conversion tests."""

import numpy as np

from projectkoios.physkit.potentials.periodic.operators import (
    PeriodicPotentialToSampledDiagonalOperatorConverter,
    PeriodicPotentialToSampledDiagonalOperatorRequest,
)
from projectkoios.physkit.potentials.periodic1d.cosine import (
    CosinePotential1D,
)
from projectkoios.physkit.units import MatrixQuantity, ScalarQuantity, Unitless


def test__action__constructs_canonical_diagonal_csr_operator() -> None:
    """Sample values occupy only the multiplication-operator diagonal."""
    potential = CosinePotential1D(
        amplitude=ScalarQuantity(2.0, Unitless()),
    )

    conversion = PeriodicPotentialToSampledDiagonalOperatorConverter().action(
        request=PeriodicPotentialToSampledDiagonalOperatorRequest(
            model=potential,
            fractional_coordinates=MatrixQuantity(
                magnitude=np.asarray([[0.0], [0.25], [0.5]]),
                unit=Unitless(),
            ),
        )
    )

    np.testing.assert_allclose(
        conversion.represented_operator.to_csr().toarray(),
        np.diag([2.0, 0.0, -2.0]),
        atol=1.0e-14,
    )
    # cos(pi/2) is a binary64 near-zero rather than an exact zero.
    assert conversion.represented_operator.nonzero_count == 3
