"""Tests for the physical closed-loop frame-transport action."""

import numpy as np
import pytest

from projectkoios.numerics.finite_periodic.meshes import CenteredUniformMesh
from projectkoios.optimization.alignment.closed_loop_polar_transport import (
    ClosedLoopPolarFrameTransportBackendFailure,
    ClosedLoopPolarFrameTransportBackendFailureReason,
    ClosedLoopPolarFrameTransportBackendOperation,
    ClosedLoopPolarFrameTransportLinkStatus,
    ClosedLoopPolarFrameTransportPrincipalPhaseBranch,
    ClosedLoopPolarFrameTransportStatus,
)
from projectkoios.physkit.solidstate.electronic import frame_transport as module
from projectkoios.physkit.solidstate.electronic.frame_transport import (
    ClosedLoopFrameTransportOutcome,
    ClosedLoopPolarFrameTransporter,
    EndpointSewingAction,
    FinitePeriodicMeshFrameInventory,
    FrameRotationAction,
    OrderedReducedPointFrameInventory,
)


def _scalar_inventory() -> OrderedReducedPointFrameInventory:
    phases = np.asarray([0.2, -1.1, 2.4], dtype=np.float64)
    frames = np.zeros((phases.size, 2, 1), dtype=np.complex128)
    frames[:, 0, 0] = np.exp(1j * phases)
    sewing = np.diag(np.asarray([np.exp(0.8j), 1.0], dtype=np.complex128))
    return OrderedReducedPointFrameInventory.from_arrays(
        reduced_points=((-0.5,), (0.0,), (0.25,)),
        frames=frames,
        endpoint_sewing_map=sewing,
        ambient_basis_identifier="ambient-plane-wave-order",
        retained_frame_identifier="isolated-band-7",
        endpoint_sewing_identifier="reciprocal-lattice-sewing",
        principal_phase_branch=(
            ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
        ),
    )


def test__transport__qualifies_scalar_route_and_reattaches_all_identities() -> None:
    source = _scalar_inventory()

    outcome = ClosedLoopPolarFrameTransporter().transport(source)

    assert isinstance(outcome, ClosedLoopFrameTransportOutcome)
    assert outcome.source is source
    assert outcome.status is ClosedLoopPolarFrameTransportStatus.QUALIFIED
    assert outcome.ambient_basis_identifier == "ambient-plane-wave-order"
    assert outcome.retained_frame_identifier == "isolated-band-7"
    assert outcome.endpoint_sewing_identifier == "reciprocal-lattice-sewing"
    assert outcome.frame_rotation_action is FrameRotationAction.RIGHT
    assert outcome.endpoint_sewing_action is EndpointSewingAction.LEFT_AMBIENT
    assert outcome.numerical_definition is outcome.numerical_result.definition
    assert outcome.numerical_definition.frames is source.frames
    assert (
        outcome.numerical_definition.endpoint_sewing_map is source.endpoint_sewing_map
    )
    assert outcome.point_count == 3
    assert outcome.ambient_dimension == 2
    assert outcome.rank == 1
    assert outcome.tolerance is source.tolerance
    assert (
        outcome.principal_phase_branch
        is ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
    )
    assert outcome.ordered_closure_endpoint_mismatch_phases == pytest.approx((-0.8,))
    assert outcome.transported_frames is outcome.numerical_result.transported_frames
    assert outcome.frame_right_actions is outcome.numerical_result.frame_right_actions
    assert outcome.frame_right_actions is not None
    np.testing.assert_allclose(
        outcome.transported_frames,
        source.frames @ outcome.frame_right_actions,
        atol=1.0e-14,
    )
    assert (
        outcome.closure_alignment_right_action
        is outcome.numerical_result.closure_alignment_right_action
    )
    assert outcome.closure_alignment_right_action is not None
    assert outcome.closure_endpoint_mismatch is not None
    np.testing.assert_allclose(
        outcome.closure_endpoint_mismatch,
        outcome.closure_alignment_right_action.conjugate().T,
        atol=1.0e-14,
    )
    assert outcome.distributed_closure_correction_step is not None
    np.testing.assert_allclose(
        np.linalg.matrix_power(
            outcome.distributed_closure_correction_step,
            outcome.point_count,
        ),
        outcome.closure_alignment_right_action,
        atol=1.0e-14,
    )
    assert outcome.closure_schur_vectors is not None
    assert outcome.maximum_absolute_entry_closure_alignment_residual is not None
    assert outcome.maximum_absolute_entry_closure_alignment_unitarity_defect is not None
    assert outcome.deficient_link_indices == ()
    assert outcome.unattempted_link_indices == ()
    assert all(
        diagnostic.status is ClosedLoopPolarFrameTransportLinkStatus.QUALIFIED
        for diagnostic in outcome.link_diagnostics
    )
    for diagnostic in outcome.link_diagnostics:
        assert diagnostic.singular_values == pytest.approx((1.0,))
        assert diagnostic.applied_rank_cutoff is not None
        assert diagnostic.rank_margin is not None
        assert diagnostic.numerical_rank == 1
        assert diagnostic.maximum_absolute_entry_unitarity_defect is not None
    assert outcome.svd_backend_identifier == "numpy.linalg.svd"
    assert outcome.svd_backend_version
    assert outcome.schur_backend_identifier == "scipy.linalg.schur"
    assert outcome.schur_backend_version
    with pytest.raises(ValueError, match="read-only"):
        outcome.transported_frames[0, 0, 0] = 0.0  # type: ignore[index]


def test__transport__qualifies_composite_finite_periodic_route() -> None:
    mesh = CenteredUniformMesh(1, 3, (0.0,))
    frames = np.broadcast_to(
        np.eye(2, dtype=np.complex128),
        (mesh.point_count, 2, 2),
    ).copy()
    mixing = np.asarray(
        [[1.0, 1.0], [1.0j, -1.0j]],
        dtype=np.complex128,
    ) / np.sqrt(2.0)
    phases = np.asarray([-0.7, 1.1], dtype=np.float64)
    sewing = mixing @ np.diag(np.exp(1j * phases)) @ mixing.conjugate().T
    source = FinitePeriodicMeshFrameInventory.from_arrays(
        mesh=mesh,
        frames=frames,
        endpoint_sewing_map=sewing,
        ambient_basis_identifier="two-orbital-order",
        retained_frame_identifier="two-band-composite",
        endpoint_sewing_identifier="mesh-period-sewing",
        principal_phase_branch=(
            ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
        ),
    )

    outcome = ClosedLoopPolarFrameTransporter().transport(source)

    assert outcome.source.mesh is mesh
    assert outcome.status is ClosedLoopPolarFrameTransportStatus.QUALIFIED
    assert outcome.closure_alignment_right_action is not None
    assert outcome.closure_endpoint_mismatch is not None
    assert outcome.distributed_closure_correction_step is not None
    np.testing.assert_allclose(
        outcome.closure_alignment_right_action,
        sewing,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        outcome.closure_endpoint_mismatch,
        sewing.conjugate().T,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        outcome.ordered_closure_endpoint_mismatch_phases,
        np.sort(-phases),
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        np.linalg.matrix_power(
            outcome.distributed_closure_correction_step,
            mesh.point_count,
        ),
        outcome.closure_alignment_right_action,
        atol=1.0e-14,
    )


def test__transport__negative_pi_open_represents_negative_one_as_positive_pi() -> None:
    frames = np.asarray([[[1.0], [0.0]]], dtype=np.complex128)
    source = OrderedReducedPointFrameInventory.from_arrays(
        reduced_points=((0.0,),),
        frames=frames,
        endpoint_sewing_map=np.diag(np.asarray([-1.0, 1.0], dtype=np.complex128)),
        ambient_basis_identifier="basis",
        retained_frame_identifier="band",
        endpoint_sewing_identifier="sewing",
        principal_phase_branch=(
            ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
        ),
    )

    outcome = ClosedLoopPolarFrameTransporter().transport(source)

    assert outcome.ordered_closure_endpoint_mismatch_phases == (np.pi,)


def test__transport__supplies_explicit_numpy_procrustes_collaborator(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    owner_solver = module.ScipyClosedLoopPolarFrameTransportSolver
    collaborators: list[object] = []

    class CapturingSolver:
        def __init__(self, procrustes_solver: object) -> None:
            collaborators.append(procrustes_solver)
            self._delegate = owner_solver(procrustes_solver)  # type: ignore[arg-type]

        def solve(self, definition: object) -> object:
            return self._delegate.solve(definition)  # type: ignore[arg-type]

    monkeypatch.setattr(
        module,
        "ScipyClosedLoopPolarFrameTransportSolver",
        CapturingSolver,
    )

    ClosedLoopPolarFrameTransporter().transport(_scalar_inventory())

    assert len(collaborators) == 1
    assert type(collaborators[0]) is module.NumpyComplexUnitaryProcrustesSolver


def test__transport__preserves_deficient_and_unattempted_links() -> None:
    frames = np.asarray(
        [
            [[1.0], [0.0]],
            [[0.0], [1.0]],
            [[0.0], [1.0]],
        ],
        dtype=np.complex128,
    )
    source = OrderedReducedPointFrameInventory.from_arrays(
        reduced_points=((0.0,), (0.25,), (0.5,)),
        frames=frames,
        endpoint_sewing_map=np.eye(2, dtype=np.complex128),
        ambient_basis_identifier="basis",
        retained_frame_identifier="band",
        endpoint_sewing_identifier="sewing",
        principal_phase_branch=(
            ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
        ),
    )

    outcome = ClosedLoopPolarFrameTransporter().transport(source)

    assert outcome.source is source
    assert outcome.status is ClosedLoopPolarFrameTransportStatus.RANK_DEFICIENT
    assert outcome.deficient_link_indices == (0,)
    assert outcome.unattempted_link_indices == (1, 2)
    assert outcome.link_diagnostics[0].singular_values == (0.0,)
    assert outcome.link_diagnostics[0].maximum_absolute_entry_unitarity_defect is None
    assert all(
        diagnostic.status is ClosedLoopPolarFrameTransportLinkStatus.UNATTEMPTED
        for diagnostic in outcome.link_diagnostics[1:]
    )
    assert outcome.transported_frames is None
    assert outcome.closure_alignment_right_action is None
    assert outcome.closure_endpoint_mismatch is None
    assert outcome.distributed_closure_correction_step is None
    assert outcome.frame_right_actions is None
    assert outcome.closure_schur_vectors is None
    assert outcome.ordered_closure_endpoint_mismatch_phases is None


def test__transport__propagates_typed_backend_failure_without_translation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    failure = ClosedLoopPolarFrameTransportBackendFailure(
        ClosedLoopPolarFrameTransportBackendOperation.PROCRUSTES_SVD,
        ClosedLoopPolarFrameTransportBackendFailureReason.NONFINITE_OUTPUT,
        1,
    )

    class FailingSolver:
        def __init__(self, _procrustes_solver: object) -> None:
            pass

        def solve(self, _definition: object) -> object:
            raise failure

    monkeypatch.setattr(
        module,
        "ScipyClosedLoopPolarFrameTransportSolver",
        FailingSolver,
    )

    with pytest.raises(ClosedLoopPolarFrameTransportBackendFailure) as captured:
        ClosedLoopPolarFrameTransporter().transport(_scalar_inventory())

    assert captured.value is failure
    assert (
        captured.value.operation
        is ClosedLoopPolarFrameTransportBackendOperation.PROCRUSTES_SVD
    )
    assert (
        captured.value.reason
        is ClosedLoopPolarFrameTransportBackendFailureReason.NONFINITE_OUTPUT
    )
    assert captured.value.link_index == 1
