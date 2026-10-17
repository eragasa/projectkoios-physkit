"""Contract tests for the unit-volume three-dimensional Bain path."""

from __future__ import annotations

import math

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice import (
    DirectLattice3D,
    UnitVolumeBainPath3D,
    UnitVolumeBainPathPoint3D,
)


def test__evaluate__retains_fcc_and_bcc_endpoint_conventions() -> None:
    """The path interpolates the represented BCT axial ratio between endpoints."""
    path = UnitVolumeBainPath3D()

    fcc = path.evaluate(0.0)
    bcc = path.evaluate(1.0)

    assert type(fcc) is UnitVolumeBainPathPoint3D
    assert type(bcc) is UnitVolumeBainPathPoint3D
    assert fcc.path_parameter == 0.0
    assert bcc.path_parameter == 1.0
    assert fcc.axial_ratio == math.sqrt(2.0)
    assert bcc.axial_ratio == 1.0
    assert fcc.longitudinal_scale / fcc.transverse_scale == math.sqrt(2.0)
    assert bcc.longitudinal_scale / bcc.transverse_scale == 1.0
    assert type(fcc.direct_lattice) is DirectLattice3D
    assert type(bcc.direct_lattice) is DirectLattice3D
    assert fcc.direct_lattice.A.flags.writeable is False
    assert bcc.direct_lattice.A.flags.writeable is False


def test__evaluate__preserves_unit_primitive_volume_on_dense_scan() -> None:
    """The represented primitive determinant remains one over the complete path."""
    path = UnitVolumeBainPath3D()
    points = tuple(path.evaluate(float(index) / 200.0) for index in range(201))
    volumes = np.asarray([point.primitive_volume for point in points])

    assert np.ptp(volumes) < 1.0e-15
    np.testing.assert_allclose(volumes, np.ones(201), rtol=0.0, atol=5.0e-16)
    for point in points:
        assert np.linalg.det(point.direct_metric.metric_tensor) == pytest.approx(
            point.primitive_volume**2,
            rel=0.0,
            abs=3.0e-15,
        )


def test__evaluate__constructs_independently_expected_primitive_basis() -> None:
    """An interior point uses the declared BCT primitive-vector convention."""
    point = UnitVolumeBainPath3D().evaluate(0.25)
    transverse = (2.0 / point.axial_ratio) ** (1.0 / 3.0)
    longitudinal = point.axial_ratio * transverse

    np.testing.assert_array_equal(
        point.direct_lattice.A,
        np.asarray(
            [
                [transverse / 2.0, -transverse / 2.0, transverse / 2.0],
                [transverse / 2.0, transverse / 2.0, -transverse / 2.0],
                [longitudinal / 2.0, longitudinal / 2.0, longitudinal / 2.0],
            ]
        ),
    )


@pytest.mark.parametrize("value", (True, 0, np.float64(0.5)))
def test__evaluate__rejects_non_builtin_float(value: object) -> None:
    """The public path parameter has one unambiguous runtime representation."""
    with pytest.raises(TypeError, match="built-in float"):
        UnitVolumeBainPath3D().evaluate(value)  # type: ignore[arg-type]


@pytest.mark.parametrize("value", (-0.01, 1.01, math.nan, math.inf))
def test__evaluate__rejects_nonfinite_or_out_of_interval_value(value: float) -> None:
    """The closed path interval excludes nonfinite and extrapolated values."""
    with pytest.raises(ValueError, match=r"finite and lie in \[0, 1\]"):
        UnitVolumeBainPath3D().evaluate(value)
