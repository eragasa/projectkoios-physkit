"""Contract tests for correlated unit-area rhombic path points."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice import (
    DirectLattice2D,
    UnitAreaRhombicPath2D,
    UnitAreaRhombicPathPoint2D,
)


def test__init__owns_lattice_and_retains_exact_path_correlation() -> None:
    """A path point owns its exact lattice rather than the caller's object."""
    source = UnitAreaRhombicPath2D().evaluate(0.5)

    point = UnitAreaRhombicPathPoint2D(
        path_parameter=source.path_parameter,
        interior_angle_radians=source.interior_angle_radians,
        basis_scale=source.basis_scale,
        direct_lattice=source.direct_lattice,
    )

    assert point.direct_lattice is not source.direct_lattice
    np.testing.assert_array_equal(point.direct_lattice.A, source.direct_lattice.A)
    assert point.direct_metric.direct_lattice is not point.direct_lattice
    np.testing.assert_array_equal(
        point.direct_metric.primitive_basis,
        point.direct_lattice.A,
    )


def test__init__rejects_inconsistent_parameter_or_lattice_state() -> None:
    """Every scalar and primitive vector must describe the same rhombic point."""
    source = UnitAreaRhombicPath2D().evaluate(0.5)
    unrelated = DirectLattice2D(
        a1=np.asarray([1.0, 0.0]),
        a2=np.asarray([0.0, 1.0]),
    )

    with pytest.raises(ValueError, match="interior_angle_radians must match"):
        UnitAreaRhombicPathPoint2D(
            path_parameter=source.path_parameter,
            interior_angle_radians=source.interior_angle_radians + 0.1,
            basis_scale=source.basis_scale,
            direct_lattice=source.direct_lattice,
        )
    with pytest.raises(ValueError, match="direct_lattice must match"):
        UnitAreaRhombicPathPoint2D(
            path_parameter=source.path_parameter,
            interior_angle_radians=source.interior_angle_radians,
            basis_scale=source.basis_scale,
            direct_lattice=unrelated,
        )
