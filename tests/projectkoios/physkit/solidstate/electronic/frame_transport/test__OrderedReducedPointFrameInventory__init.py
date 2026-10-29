"""Tests for ordered reduced-point frame inventories."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from projectkoios.physkit.solidstate.electronic.frame_transport import (
    ClosedLoopPolarFrameTransportPrincipalPhaseBranch,
    EndpointSewingAction,
    FrameRotationAction,
    OrderedReducedPointFrameInventory,
)


def _frames() -> np.ndarray:
    return np.asarray(
        [
            [[1.0], [0.0]],
            [[np.exp(0.2j)], [0.0]],
        ],
        dtype=np.complex128,
    )


def test__init__retains_exact_route_identities_and_immutable_storage() -> None:
    frames = _frames()
    sewing = np.eye(2, dtype=np.complex128)
    points = ((-0.5,), (0.0,))

    inventory = OrderedReducedPointFrameInventory.from_arrays(
        reduced_points=points,
        frames=frames,
        endpoint_sewing_map=sewing,
        ambient_basis_identifier="plane-wave-order",
        retained_frame_identifier="bands-4-through-4",
        endpoint_sewing_identifier="reciprocal-period-sewing",
        principal_phase_branch=(
            ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
        ),
    )
    frames[0, 0, 0] = 0.0
    sewing[0, 0] = 0.0

    assert inventory.reduced_points is points
    assert inventory.frames[0, 0, 0] == 1.0 + 0.0j
    assert inventory.endpoint_sewing_map[0, 0] == 1.0 + 0.0j
    assert inventory.frame_rotation_action is FrameRotationAction.RIGHT
    assert inventory.endpoint_sewing_action is EndpointSewingAction.LEFT_AMBIENT
    assert inventory.ambient_basis_identifier == "plane-wave-order"
    assert inventory.retained_frame_identifier == "bands-4-through-4"
    assert inventory.endpoint_sewing_identifier == "reciprocal-period-sewing"
    assert (
        inventory.principal_phase_branch
        is ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
    )
    with pytest.raises(ValueError, match="read-only"):
        inventory.frames[0, 0, 0] = 0.0
    with pytest.raises(FrozenInstanceError):
        inventory.retained_frame_identifier = "other"  # type: ignore[misc]


def test__init__rejects_route_frame_mismatch_without_route_policy() -> None:
    with pytest.raises(ValueError, match="length must equal"):
        OrderedReducedPointFrameInventory.from_arrays(
            reduced_points=((0.0,),),
            frames=_frames(),
            endpoint_sewing_map=np.eye(2, dtype=np.complex128),
            ambient_basis_identifier="basis",
            retained_frame_identifier="frame",
            endpoint_sewing_identifier="sewing",
            principal_phase_branch=(
                ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
            ),
        )

    with pytest.raises(ValueError, match="same dimension"):
        OrderedReducedPointFrameInventory.from_arrays(
            reduced_points=((0.0,), (0.0, 0.5)),
            frames=_frames(),
            endpoint_sewing_map=np.eye(2, dtype=np.complex128),
            ambient_basis_identifier="basis",
            retained_frame_identifier="frame",
            endpoint_sewing_identifier="sewing",
            principal_phase_branch=(
                ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open
            ),
        )
