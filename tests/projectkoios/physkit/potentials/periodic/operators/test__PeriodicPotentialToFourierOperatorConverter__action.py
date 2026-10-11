"""Reciprocal-mode periodic-potential operator conversion tests."""

import numpy as np

from projectkoios.physkit.potentials.periodic.operators import (
    PeriodicPotentialToFourierOperatorConverter,
    PeriodicPotentialToFourierOperatorRequest,
)
from projectkoios.physkit.potentials.periodic1d.cosine import (
    CosinePotential1D,
)
from projectkoios.physkit.units import ScalarQuantity, Unitless


def test__action__constructs_hermitian_mode_coupling_operator() -> None:
    """A real cosine couples reciprocal modes separated by one harmonic."""
    potential = CosinePotential1D(
        amplitude=ScalarQuantity(4.0, Unitless()),
    )

    conversion = PeriodicPotentialToFourierOperatorConverter().action(
        request=PeriodicPotentialToFourierOperatorRequest(
            model=potential,
            mode_indices=((-1,), (0,), (1,)),
        )
    )

    represented = conversion.represented_operator.to_csr().toarray()
    np.testing.assert_allclose(
        represented,
        np.asarray(
            [
                [0.0, 2.0, 0.0],
                [2.0, 0.0, 2.0],
                [0.0, 2.0, 0.0],
            ]
        ),
    )
    np.testing.assert_allclose(represented, represented.conj().T)
    assert conversion.represented_operator.nonzero_count == 4
