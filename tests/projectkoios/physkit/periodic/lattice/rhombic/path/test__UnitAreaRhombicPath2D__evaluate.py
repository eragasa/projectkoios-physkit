"""Contract tests for the unit-area two-dimensional rhombic path."""

from __future__ import annotations

import math

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice import (
    DirectLattice2D,
    UnitAreaRhombicPath2D,
    UnitAreaRhombicPathPoint2D,
)


def test__evaluate__retains_triangular_and_square_endpoint_conventions() -> None:
    """The path interpolates the primitive angle between declared endpoints."""
    path = UnitAreaRhombicPath2D()

    triangular = path.evaluate(0.0)
    square = path.evaluate(1.0)

    assert type(triangular) is UnitAreaRhombicPathPoint2D
    assert type(square) is UnitAreaRhombicPathPoint2D
    assert triangular.interior_angle_radians == math.pi / 3.0
    assert square.interior_angle_radians == math.pi / 2.0
    assert type(triangular.direct_lattice) is DirectLattice2D
    assert type(square.direct_lattice) is DirectLattice2D
    np.testing.assert_allclose(
        square.direct_metric.metric_tensor, np.eye(2), atol=2e-16
    )


def test__evaluate__preserves_unit_area_on_dense_scan() -> None:
    """The represented primitive determinant remains one over the complete path."""
    points = tuple(
        UnitAreaRhombicPath2D().evaluate(float(index) / 200.0) for index in range(201)
    )
    areas = np.asarray([point.primitive_area for point in points])

    assert np.ptp(areas) < 6.0e-16
    np.testing.assert_allclose(areas, np.ones(201), rtol=0.0, atol=4.0e-16)
    for point in points:
        angle = point.interior_angle_radians
        expected_metric = np.asarray(
            ((1.0, math.cos(angle)), (math.cos(angle), 1.0))
        ) / math.sin(angle)
        np.testing.assert_allclose(
            point.direct_metric.metric_tensor,
            expected_metric,
            rtol=0.0,
            atol=5.0e-16,
        )


def test__evaluate__constructs_declared_canonical_direct_basis() -> None:
    """An interior point uses first-vector-axis and positive-angle orientation."""
    point = UnitAreaRhombicPath2D().evaluate(0.4)
    angle = point.interior_angle_radians
    scale = point.basis_scale

    np.testing.assert_array_equal(
        point.direct_lattice.A,
        [[scale, scale * math.cos(angle)], [0.0, scale * math.sin(angle)]],
    )


@pytest.mark.parametrize("value", (True, 0, np.float64(0.5)))
def test__evaluate__rejects_non_builtin_float(value: object) -> None:
    """The public path parameter has one unambiguous runtime representation."""
    with pytest.raises(TypeError, match="built-in float"):
        UnitAreaRhombicPath2D().evaluate(value)  # type: ignore[arg-type]


@pytest.mark.parametrize("value", (-0.01, 1.01, math.nan, math.inf))
def test__evaluate__rejects_nonfinite_or_out_of_interval_value(value: float) -> None:
    """The closed path interval excludes nonfinite and extrapolated values."""
    with pytest.raises(ValueError, match=r"finite and lie in \[0, 1\]"):
        UnitAreaRhombicPath2D().evaluate(value)
