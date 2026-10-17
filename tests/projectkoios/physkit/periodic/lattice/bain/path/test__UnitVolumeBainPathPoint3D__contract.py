"""Contract tests for correlated unit-volume Bain path points."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice import (
    DirectLattice3D,
    UnitVolumeBainPath3D,
    UnitVolumeBainPathPoint3D,
)


def test__init__owns_lattice_and_retains_exact_path_correlation() -> None:
    """A path point owns its exact lattice rather than the caller's object."""
    source = UnitVolumeBainPath3D().evaluate(0.5)

    point = UnitVolumeBainPathPoint3D(
        path_parameter=source.path_parameter,
        axial_ratio=source.axial_ratio,
        transverse_scale=source.transverse_scale,
        longitudinal_scale=source.longitudinal_scale,
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
    """Every scalar and primitive vector must describe the same Bain point."""
    source = UnitVolumeBainPath3D().evaluate(0.5)
    unrelated = DirectLattice3D(
        a1=np.asarray([1.0, 0.0, 0.0]),
        a2=np.asarray([0.0, 1.0, 0.0]),
        a3=np.asarray([0.0, 0.0, 1.0]),
    )

    with pytest.raises(ValueError, match="axial_ratio must match"):
        UnitVolumeBainPathPoint3D(
            path_parameter=source.path_parameter,
            axial_ratio=source.axial_ratio + 0.1,
            transverse_scale=source.transverse_scale,
            longitudinal_scale=source.longitudinal_scale,
            direct_lattice=source.direct_lattice,
        )
    with pytest.raises(ValueError, match="direct_lattice must match"):
        UnitVolumeBainPathPoint3D(
            path_parameter=source.path_parameter,
            axial_ratio=source.axial_ratio,
            transverse_scale=source.transverse_scale,
            longitudinal_scale=source.longitudinal_scale,
            direct_lattice=unrelated,
        )
