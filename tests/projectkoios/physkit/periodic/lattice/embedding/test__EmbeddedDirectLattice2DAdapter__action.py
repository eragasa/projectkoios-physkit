"""Numerical verification for ambient-to-intrinsic lattice adaptation."""

from __future__ import annotations

import numpy as np

from projectkoios.physkit.periodic.finite_difference.grid import (
    PeriodicFiniteDifferenceGrid,
)
from projectkoios.physkit.periodic.finite_difference.laplacian import (
    BlochPeriodicLaplacianConstructionRequest,
    BlochPeriodicLaplacianConstructor,
)
from projectkoios.physkit.periodic.lattice.boundary_phases import BoundaryTwistLift
from projectkoios.physkit.periodic.lattice.embedding import (
    EmbeddedDirectLattice2D,
    EmbeddedDirectLattice2DAdaptation,
    EmbeddedDirectLattice2DAdaptationRequest,
    EmbeddedDirectLattice2DAdapter,
)
from projectkoios.physkit.periodic.lattice.finite_domain import (
    FinitePeriodicDomain,
    LatticeDimension,
)
from projectkoios.physkit.periodic.lattice.lattice2d import DirectLattice2D
from projectkoios.physkit.periodic.lattice.metric import DirectLatticeMetric
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)


def _adapt(
    basis: np.ndarray,
    unit: PhysicalUnit | Unitless | None = None,
) -> EmbeddedDirectLattice2DAdaptation:
    selected_unit = Unitless() if unit is None else unit
    return EmbeddedDirectLattice2DAdapter().action(
        request=EmbeddedDirectLattice2DAdaptationRequest(
            embedded_lattice=EmbeddedDirectLattice2D(
                ambient_basis=MatrixQuantity(basis, selected_unit)
            )
        )
    )


def _laplacian(
    adaptation: EmbeddedDirectLattice2DAdaptation,
    *,
    shape: tuple[int, int] = (5, 6),
    twist: tuple[float, float] = (0.17, 0.31),
):
    grid = PeriodicFiniteDifferenceGrid(
        metric=adaptation.intrinsic_metric,
        domain=FinitePeriodicDomain(LatticeDimension.TWO, shape),
        direct_basis_unit=adaptation.request.embedded_lattice.basis_unit,
    )
    return BlochPeriodicLaplacianConstructor().action(
        request=BlochPeriodicLaplacianConstructionRequest(
            grid=grid,
            twist=BoundaryTwistLift(LatticeDimension.TWO, twist),
        )
    )


def _intrinsic_laplacian(
    basis: np.ndarray,
    unit: PhysicalUnit | Unitless,
    *,
    shape: tuple[int, int] = (5, 6),
    twist: tuple[float, float] = (0.17, 0.31),
):
    lattice = DirectLattice2D(a1=basis[:, 0], a2=basis[:, 1])
    grid = PeriodicFiniteDifferenceGrid(
        metric=DirectLatticeMetric(lattice),
        domain=FinitePeriodicDomain(LatticeDimension.TWO, shape),
        direct_basis_unit=unit,
    )
    return BlochPeriodicLaplacianConstructor().action(
        request=BlochPeriodicLaplacianConstructionRequest(
            grid=grid,
            twist=BoundaryTwistLift(LatticeDimension.TWO, twist),
        )
    )


def test__action__constructs_metric_dual_frame_and_projector_identities() -> None:
    basis = np.array(
        [
            [2.0, 0.4],
            [0.3, 1.7],
            [1.1, -0.2],
        ]
    )
    result = _adapt(basis)
    metric = basis.T @ basis
    frame = result.orthonormal_frame.magnitude
    projector = result.ambient_projector.magnitude

    np.testing.assert_allclose(result.induced_metric.magnitude, metric, rtol=2e-15)
    np.testing.assert_allclose(
        result.inverse_induced_metric.magnitude,
        np.linalg.inv(metric),
        rtol=3e-15,
    )
    np.testing.assert_allclose(frame.T @ frame, np.identity(2), atol=3e-15)
    np.testing.assert_allclose(projector, frame @ frame.T, atol=3e-15)
    np.testing.assert_allclose(projector @ projector, projector, atol=3e-15)
    np.testing.assert_allclose(
        basis.T @ result.dual_reciprocal_basis.magnitude,
        2.0 * np.pi * np.identity(2),
        rtol=3e-15,
        atol=3e-15,
    )
    np.testing.assert_allclose(
        result.intrinsic_basis.magnitude.T @ result.intrinsic_basis.magnitude,
        metric,
        rtol=3e-15,
        atol=3e-15,
    )


def test__action__maps_fractional_intrinsic_and_ambient_vectors() -> None:
    basis = np.array(
        [
            [2.0, 0.4],
            [0.3, 1.7],
            [1.1, -0.2],
        ]
    )
    result = _adapt(basis, PhysicalUnit("angstrom"))
    fractional = np.array([[0.2, -0.3], [1.1, 0.4]])
    ambient = result.fractional_to_ambient(fractional)

    np.testing.assert_allclose(ambient, fractional @ basis.T, rtol=2e-15)
    np.testing.assert_allclose(
        result.project_ambient_to_fractional(ambient),
        fractional,
        rtol=3e-15,
        atol=3e-15,
    )

    intrinsic = fractional @ result.intrinsic_basis.magnitude.T
    np.testing.assert_allclose(
        result.intrinsic_to_ambient(intrinsic),
        ambient,
        rtol=3e-15,
        atol=3e-15,
    )
    np.testing.assert_allclose(
        result.project_ambient_to_intrinsic(ambient),
        intrinsic,
        rtol=3e-15,
        atol=3e-15,
    )

    normal_component = np.array([0.3, -0.5, 0.8])
    normal_component -= result.ambient_projector.magnitude @ normal_component
    displaced = ambient + normal_component
    np.testing.assert_allclose(
        result.project_ambient_to_fractional(displaced),
        fractional,
        rtol=3e-15,
        atol=3e-15,
    )


def test__action__is_invariant_under_ambient_rotation_and_reflection() -> None:
    basis = np.array(
        [
            [2.0, 0.4],
            [0.3, 1.7],
            [1.1, -0.2],
        ]
    )
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    reflection = np.diag([-1.0, 1.0, 1.0])
    original = _adapt(basis)
    rotated = _adapt(rotation @ basis)
    reflected = _adapt(reflection @ basis)

    for transformed, ambient_operation in (
        (rotated, rotation),
        (reflected, reflection),
    ):
        np.testing.assert_allclose(
            transformed.induced_metric.magnitude,
            original.induced_metric.magnitude,
            rtol=3e-15,
            atol=3e-15,
        )
        np.testing.assert_allclose(
            transformed.intrinsic_basis.magnitude,
            original.intrinsic_basis.magnitude,
            rtol=3e-15,
            atol=3e-15,
        )
        np.testing.assert_allclose(
            transformed.orthonormal_frame.magnitude,
            ambient_operation @ original.orthonormal_frame.magnitude,
            rtol=3e-15,
            atol=3e-15,
        )
        np.testing.assert_allclose(
            transformed.ambient_projector.magnitude,
            ambient_operation
            @ original.ambient_projector.magnitude
            @ ambient_operation.T,
            rtol=3e-15,
            atol=3e-15,
        )
        np.testing.assert_allclose(
            _laplacian(transformed).represented_laplacian.to_csr().toarray(),
            _laplacian(original).represented_laplacian.to_csr().toarray(),
            rtol=2e-14,
            atol=2e-14,
        )


def test__action__matches_explicit_intrinsic_route_with_bloch_seams() -> None:
    basis = np.array(
        [
            [2.0, 0.4],
            [0.3, 1.7],
            [1.1, -0.2],
        ]
    )
    embedded = _adapt(basis)
    intrinsic_basis = embedded.intrinsic_basis.magnitude
    direct_intrinsic = _intrinsic_laplacian(intrinsic_basis, Unitless())

    np.testing.assert_allclose(
        intrinsic_basis.T @ intrinsic_basis,
        embedded.induced_metric.magnitude,
        rtol=3e-15,
        atol=3e-15,
    )
    np.testing.assert_allclose(
        direct_intrinsic.represented_laplacian.to_csr().toarray(),
        _laplacian(embedded).represented_laplacian.to_csr().toarray(),
        rtol=2e-14,
        atol=2e-14,
    )


def test__action__preserves_physical_scale_and_unitless_support() -> None:
    basis_angstrom = np.array(
        [
            [3.0, 0.8],
            [0.2, 2.1],
            [1.0, -0.4],
        ]
    )
    angstrom = _adapt(basis_angstrom, PhysicalUnit("angstrom"))
    nanometer = _adapt(basis_angstrom / 10.0, PhysicalUnit("nanometer"))
    unitless = _adapt(basis_angstrom, Unitless())

    basis_nm_to_angstrom = MODEL_SYSTEM_UNIT_CONVERTER.convert_matrix(
        nanometer.request.embedded_lattice.ambient_basis,
        PhysicalUnit("angstrom"),
    )
    metric_nm_to_angstrom = MODEL_SYSTEM_UNIT_CONVERTER.convert_matrix(
        nanometer.induced_metric,
        PhysicalUnit("angstrom ** 2"),
    )
    dual_nm_to_angstrom = MODEL_SYSTEM_UNIT_CONVERTER.convert_matrix(
        nanometer.dual_reciprocal_basis,
        PhysicalUnit("1 / angstrom"),
    )
    laplacian_nm_to_angstrom = (
        MODEL_SYSTEM_UNIT_CONVERTER.convert_complex_sparse_matrix(
            _laplacian(nanometer).represented_laplacian,
            PhysicalUnit("1 / angstrom ** 2"),
        )
    )

    np.testing.assert_allclose(
        basis_nm_to_angstrom.magnitude,
        angstrom.request.embedded_lattice.ambient_basis.magnitude,
        rtol=2e-15,
        atol=2e-15,
    )
    np.testing.assert_allclose(
        metric_nm_to_angstrom.magnitude,
        angstrom.induced_metric.magnitude,
        rtol=2e-14,
        atol=2e-14,
    )
    np.testing.assert_allclose(
        dual_nm_to_angstrom.magnitude,
        angstrom.dual_reciprocal_basis.magnitude,
        rtol=2e-14,
        atol=2e-14,
    )
    np.testing.assert_allclose(
        laplacian_nm_to_angstrom.to_csr().toarray(),
        _laplacian(angstrom).represented_laplacian.to_csr().toarray(),
        rtol=3e-14,
        atol=3e-14,
    )
    assert isinstance(unitless.induced_metric.unit, Unitless)
    assert isinstance(unitless.dual_reciprocal_basis.unit, Unitless)
    assert isinstance(_laplacian(unitless).represented_laplacian.unit, Unitless)
