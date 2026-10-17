"""Software and numerical verification for nonorthogonal Bloch Laplacians."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from projectkoios.physkit.periodic.finite_difference.grid import (
    PeriodicFiniteDifferenceGrid,
)
from projectkoios.physkit.periodic.finite_difference.laplacian import (
    BlochPeriodicLaplacian,
    BlochPeriodicLaplacianConstructionRequest,
    BlochPeriodicLaplacianConstructor,
)
from projectkoios.physkit.periodic.lattice.boundary_phases import BoundaryTwistLift
from projectkoios.physkit.periodic.lattice.finite_domain import (
    FinitePeriodicDomain,
    LatticeDimension,
)
from projectkoios.physkit.periodic.lattice.lattice2d import DirectLattice2D
from projectkoios.physkit.periodic.lattice.lattice3d import DirectLattice3D
from projectkoios.physkit.periodic.lattice.metric import DirectLatticeMetric
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    Unitless,
)


def _construct(
    basis: np.ndarray,
    shape: tuple[int, ...],
    twist_turns: tuple[float, ...],
    *,
    unit: PhysicalUnit | Unitless | None = None,
) -> BlochPeriodicLaplacian:
    dimension = len(shape)
    if dimension == 2:
        lattice = DirectLattice2D(a1=basis[:, 0], a2=basis[:, 1])
        lattice_dimension = LatticeDimension.TWO
        closed_shape = (shape[0], shape[1])
        closed_twist = (twist_turns[0], twist_turns[1])
    else:
        lattice = DirectLattice3D(
            a1=basis[:, 0],
            a2=basis[:, 1],
            a3=basis[:, 2],
        )
        lattice_dimension = LatticeDimension.THREE
        closed_shape = (shape[0], shape[1], shape[2])
        closed_twist = (twist_turns[0], twist_turns[1], twist_turns[2])
    grid = PeriodicFiniteDifferenceGrid(
        metric=DirectLatticeMetric(lattice),
        domain=FinitePeriodicDomain(lattice_dimension, closed_shape),
        direct_basis_unit=Unitless() if unit is None else unit,
    )
    return BlochPeriodicLaplacianConstructor().action(
        request=BlochPeriodicLaplacianConstructionRequest(
            grid=grid,
            twist=BoundaryTwistLift(lattice_dimension, closed_twist),
        )
    )


def _mode(
    shape: tuple[int, ...],
    mode_indices: tuple[int, ...],
    twist_turns: tuple[float, ...],
) -> np.ndarray:
    coordinates = np.indices(shape, dtype=np.float64)
    phase = np.zeros(shape, dtype=np.float64)
    for axis, extent in enumerate(shape):
        phase += (
            2.0
            * np.pi
            * (mode_indices[axis] + twist_turns[axis])
            * coordinates[axis]
            / extent
        )
    return np.exp(1.0j * phase).ravel(order="C")


def _discrete_symbol(
    inverse_metric: np.ndarray,
    shape: tuple[int, ...],
    mode_indices: tuple[int, ...],
    twist_turns: tuple[float, ...],
) -> float:
    theta = np.asarray(
        [
            2.0 * np.pi * (mode + twist) / extent
            for mode, twist, extent in zip(
                mode_indices,
                twist_turns,
                shape,
                strict=True,
            )
        ],
        dtype=np.float64,
    )
    spacings = 1.0 / np.asarray(shape, dtype=np.float64)
    result = 0.0
    for axis in range(len(shape)):
        result -= (
            4.0
            * inverse_metric[axis, axis]
            * np.sin(theta[axis] / 2.0) ** 2
            / spacings[axis] ** 2
        )
    for first_axis in range(len(shape)):
        for second_axis in range(first_axis + 1, len(shape)):
            result -= (
                2.0
                * inverse_metric[first_axis, second_axis]
                * np.sin(theta[first_axis])
                * np.sin(theta[second_axis])
                / (spacings[first_axis] * spacings[second_axis])
            )
    return float(result)


def _axis_second_derivative(extent: int, twist: float) -> sparse.csr_array:
    spacing = 1.0 / extent
    represented = sparse.diags(
        (
            np.ones(extent - 1),
            -2.0 * np.ones(extent),
            np.ones(extent - 1),
        ),
        offsets=(-1, 0, 1),
        format="lil",
        dtype=np.complex128,
    )
    phase = np.exp(2.0j * np.pi * twist)
    represented[0, -1] = np.conjugate(phase)
    represented[-1, 0] = phase
    return sparse.csr_array(represented) / spacing**2


def test__action__uses_c_order_and_multi_axis_quotient_seam_phases() -> None:
    basis = np.array([[2.0, 0.7], [0.0, 1.3]])
    shape = (4, 5)
    twist = (1.25, -0.2)
    result = _construct(basis, shape, twist)
    represented = result.represented_laplacian.to_csr()
    inverse_metric = result.request.grid.metric.inverse_metric_tensor
    h0, h1 = result.request.grid.fractional_spacings
    mixed_scale = inverse_metric[0, 1] / (2.0 * h0 * h1)

    source_last_last = (shape[0] - 1) * shape[1] + shape[1] - 1
    target_zero_zero = 0
    expected_positive_corner = mixed_scale * np.exp(2.0j * np.pi * (0.25 + 0.8))
    np.testing.assert_allclose(
        represented[source_last_last, target_zero_zero],
        expected_positive_corner,
        rtol=2e-15,
        atol=2e-15,
    )

    source_zero_last = shape[1] - 1
    target_last_zero = (shape[0] - 1) * shape[1]
    expected_negative_positive_corner = -mixed_scale * np.exp(
        2.0j * np.pi * (-0.25 + 0.8)
    )
    np.testing.assert_allclose(
        represented[source_zero_last, target_last_zero],
        expected_negative_positive_corner,
        rtol=2e-15,
        atol=2e-15,
    )
    assert result.twist_fiber.representative.turns == (0.25, 0.8)


def test__action__is_hermitian_and_zero_twist_has_constant_null_mode() -> None:
    basis = np.array([[2.0, 0.8], [0.0, 1.4]])
    twisted = _construct(basis, (5, 6), (0.17, 0.31))
    represented = twisted.represented_laplacian.to_csr()
    residual = represented - represented.conjugate().T

    assert residual.nnz == 0 or np.max(np.abs(residual.data)) < 2.0e-13

    periodic = _construct(basis, (5, 6), (0.0, 0.0))
    periodic_matrix = periodic.represented_laplacian.to_csr()
    np.testing.assert_allclose(
        periodic_matrix @ np.ones(periodic.request.grid.point_count),
        0.0,
        rtol=0.0,
        atol=2.0e-13,
    )


def test__action__reduces_to_orthogonal_kronecker_sum() -> None:
    basis = np.diag([2.0, 3.0])
    shape = (4, 5)
    twist = (0.13, 0.27)
    result = _construct(basis, shape, twist)
    represented = result.represented_laplacian.to_csr()
    first = _axis_second_derivative(shape[0], twist[0])
    second = _axis_second_derivative(shape[1], twist[1])
    expected = 0.25 * sparse.kron(
        first,
        sparse.identity(shape[1], dtype=np.complex128, format="csr"),
        format="csr",
    ) + (1.0 / 9.0) * sparse.kron(
        sparse.identity(shape[0], dtype=np.complex128, format="csr"),
        second,
        format="csr",
    )

    np.testing.assert_allclose(
        represented.toarray(),
        expected.toarray(),
        rtol=2e-15,
        atol=2e-15,
    )


def test__action__matches_skew_two_dimensional_plane_wave_symbol() -> None:
    basis = np.array([[2.1, 0.9], [0.2, 1.7]])
    shape = (7, 8)
    twist = (0.19, 0.36)
    mode_indices = (1, -2)
    result = _construct(basis, shape, twist)
    represented = result.represented_laplacian.to_csr()
    eigenvector = _mode(shape, mode_indices, twist)
    expected = _discrete_symbol(
        result.request.grid.metric.inverse_metric_tensor,
        shape,
        mode_indices,
        twist,
    )

    np.testing.assert_allclose(
        represented @ eigenvector,
        expected * eigenvector,
        rtol=3e-14,
        atol=3e-13,
    )


def test__action__matches_nonorthogonal_three_dimensional_plane_wave_symbol() -> None:
    basis = 0.5 * np.array(
        [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
        ]
    )
    shape = (5, 6, 7)
    twist = (0.13, 0.21, 0.34)
    mode_indices = (1, -1, 2)
    result = _construct(basis, shape, twist)
    represented = result.represented_laplacian.to_csr()
    eigenvector = _mode(shape, mode_indices, twist)
    expected = _discrete_symbol(
        result.request.grid.metric.inverse_metric_tensor,
        shape,
        mode_indices,
        twist,
    )

    assert np.any(np.triu(result.request.grid.metric.inverse_metric_tensor, k=1) != 0.0)
    np.testing.assert_allclose(
        represented @ eigenvector,
        expected * eigenvector,
        rtol=7e-14,
        atol=2e-12,
    )


def test__action__has_second_order_plane_wave_convergence() -> None:
    basis = np.array([[1.8, 0.7], [0.1, 1.4]])
    twist = (0.17, 0.29)
    mode_indices = (1, -1)
    continuum_wave_vector = (
        2.0
        * np.pi
        * (
            np.asarray(mode_indices, dtype=np.float64)
            + np.asarray(twist, dtype=np.float64)
        )
    )
    inverse_metric = np.linalg.inv(basis.T @ basis)
    continuum_symbol = -float(
        continuum_wave_vector @ inverse_metric @ continuum_wave_vector
    )
    errors: list[float] = []

    for extent in (8, 16, 32):
        shape = (extent, extent)
        result = _construct(basis, shape, twist)
        vector = _mode(shape, mode_indices, twist)
        image = result.represented_laplacian.to_csr() @ vector
        observed = float(np.vdot(vector, image).real / np.vdot(vector, vector).real)
        errors.append(abs(observed - continuum_symbol))

    assert errors[0] / errors[1] > 3.5
    assert errors[1] / errors[2] > 3.8


def test__action__preserves_physical_scale_across_length_units() -> None:
    basis_angstrom = np.array([[3.0, 1.0], [0.0, 2.0]])
    shape = (4, 5)
    twist = (0.11, 0.37)
    angstrom = _construct(
        basis_angstrom,
        shape,
        twist,
        unit=PhysicalUnit("angstrom"),
    )
    nanometer = _construct(
        basis_angstrom / 10.0,
        shape,
        twist,
        unit=PhysicalUnit("nanometer"),
    )
    converted = ModelSystemUnitConverter.convert_complex_sparse_matrix(
        nanometer.represented_laplacian,
        PhysicalUnit("1 / angstrom ** 2"),
    )

    assert angstrom.represented_laplacian.unit == PhysicalUnit("1 / (angstrom) ** 2")
    np.testing.assert_allclose(
        converted.to_csr().toarray(),
        angstrom.represented_laplacian.to_csr().toarray(),
        rtol=2e-14,
        atol=2e-14,
    )


def test__action__returns_canonical_immutable_csr_quantity() -> None:
    result = _construct(
        np.array([[2.0, 0.6], [0.0, 1.5]]),
        (4, 5),
        (0.2, 0.3),
    )
    quantity = result.represented_laplacian
    original_data = quantity.data.copy()
    mutable_copy = quantity.to_csr()
    mutable_copy.data[:] = 0.0

    assert quantity.shape == (20, 20)
    assert quantity.data.flags.writeable is False
    assert quantity.column_indices.flags.writeable is False
    assert quantity.row_offsets.flags.writeable is False
    assert np.array_equal(quantity.data, original_data)
    for row in range(quantity.shape[0]):
        columns = quantity.column_indices[
            quantity.row_offsets[row] : quantity.row_offsets[row + 1]
        ]
        assert np.all(columns[1:] > columns[:-1])
