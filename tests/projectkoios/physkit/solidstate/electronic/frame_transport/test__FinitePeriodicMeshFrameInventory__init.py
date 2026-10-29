"""Tests for finite-periodic mesh frame inventories."""

import numpy as np
import pytest

from projectkoios.numerics.finite_periodic.meshes import CenteredUniformMesh
from projectkoios.physkit.solidstate.electronic.frame_transport import (
    ClosedLoopPolarFrameTransportPrincipalPhaseBranch,
    FinitePeriodicMeshFrameInventory,
)


def test__init__retains_exact_mesh_identity() -> None:
    mesh = CenteredUniformMesh(1, 3, (0.25,))
    frames = np.broadcast_to(
        np.eye(2, dtype=np.complex128),
        (mesh.point_count, 2, 2),
    ).copy()

    inventory = FinitePeriodicMeshFrameInventory.from_arrays(
        mesh=mesh,
        frames=frames,
        endpoint_sewing_map=np.eye(2, dtype=np.complex128),
        ambient_basis_identifier="ordered-orbitals",
        retained_frame_identifier="two-band-composite",
        endpoint_sewing_identifier="periodic-sewing",
        principal_phase_branch=(
            ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
        ),
    )

    assert inventory.mesh is mesh
    assert inventory.point_count == mesh.point_count
    assert inventory.ambient_dimension == 2
    assert inventory.rank == 2


@pytest.mark.parametrize(
    "mesh",
    [
        CenteredUniformMesh(2, 2, (0.0, 0.0)),
        CenteredUniformMesh(3, 2, (0.0, 0.0, 0.0)),
    ],
)
def test__init__rejects_mesh_that_does_not_define_one_closed_route(
    mesh: CenteredUniformMesh,
) -> None:
    frames = np.zeros((mesh.point_count, 2, 1), dtype=np.complex128)
    frames[:, 0, 0] = 1.0

    with pytest.raises(ValueError, match="dimension must equal one"):
        FinitePeriodicMeshFrameInventory.from_arrays(
            mesh=mesh,
            frames=frames,
            endpoint_sewing_map=np.eye(2, dtype=np.complex128),
            ambient_basis_identifier="basis",
            retained_frame_identifier="frame",
            endpoint_sewing_identifier="sewing",
            principal_phase_branch=(
                ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
            ),
        )


def test__init__rejects_frame_count_not_owned_by_mesh() -> None:
    mesh = CenteredUniformMesh(1, 3, (0.0,))
    frames = np.broadcast_to(
        np.eye(2, dtype=np.complex128),
        (2, 2, 2),
    ).copy()

    with pytest.raises(ValueError, match="mesh point count"):
        FinitePeriodicMeshFrameInventory.from_arrays(
            mesh=mesh,
            frames=frames,
            endpoint_sewing_map=np.eye(2, dtype=np.complex128),
            ambient_basis_identifier="basis",
            retained_frame_identifier="frame",
            endpoint_sewing_identifier="sewing",
            principal_phase_branch=(
                ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
            ),
        )
