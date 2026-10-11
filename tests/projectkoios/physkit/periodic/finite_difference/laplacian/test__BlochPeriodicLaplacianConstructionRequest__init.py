"""Contract tests for Bloch-periodic Laplacian requests."""

import numpy as np
import pytest

from projectkoios.physkit.periodic.finite_difference.grid import (
    PeriodicFiniteDifferenceGrid,
)
from projectkoios.physkit.periodic.finite_difference.laplacian import (
    BlochPeriodicLaplacianConstructionRequest,
)
from projectkoios.physkit.periodic.lattice.boundary_phases import BoundaryTwistLift
from projectkoios.physkit.periodic.lattice.finite_domain import (
    FinitePeriodicDomain,
    LatticeDimension,
)
from projectkoios.physkit.periodic.lattice.lattice2d import DirectLattice2D
from projectkoios.physkit.periodic.lattice.metric import DirectLatticeMetric
from projectkoios.physkit.units import Unitless


def test__init__rejects_twist_dimension_mismatch() -> None:
    metric = DirectLatticeMetric(
        DirectLattice2D(
            a1=np.array([1.0, 0.0]),
            a2=np.array([0.2, 1.0]),
        )
    )
    grid = PeriodicFiniteDifferenceGrid(
        metric=metric,
        domain=FinitePeriodicDomain(LatticeDimension.TWO, (3, 4)),
        direct_basis_unit=Unitless(),
    )

    with pytest.raises(ValueError, match="dimensions must agree"):
        BlochPeriodicLaplacianConstructionRequest(
            grid=grid,
            twist=BoundaryTwistLift(
                LatticeDimension.THREE,
                (0.0, 0.0, 0.0),
            ),
        )
