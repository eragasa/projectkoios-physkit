"""Public API checks for physical closed-loop frame transport."""

from projectkoios.optimization.alignment import (
    closed_loop_polar_transport as contract,
)
from projectkoios.physkit.solidstate import electronic
from projectkoios.physkit.solidstate.electronic import frame_transport


def test__public_api__exports_wrapper_and_exact_owner_evidence_types() -> None:
    wrapper_names = (
        "ClosedLoopFrameTransportOutcome",
        "ClosedLoopPolarFrameTransporter",
        "EndpointSewingAction",
        "FinitePeriodicMeshFrameInventory",
        "FrameRotationAction",
        "OrderedReducedPointFrameInventory",
    )
    owner_names = (
        "ClosedLoopPolarFrameTransportBackendFailure",
        "ClosedLoopPolarFrameTransportBackendFailureReason",
        "ClosedLoopPolarFrameTransportBackendOperation",
        "ClosedLoopPolarFrameTransportLinkDiagnostic",
        "ClosedLoopPolarFrameTransportLinkKind",
        "ClosedLoopPolarFrameTransportLinkStatus",
        "ClosedLoopPolarFrameTransportPrincipalPhaseBranch",
        "ClosedLoopPolarFrameTransportStatus",
    )

    for name in wrapper_names:
        assert getattr(electronic, name) is getattr(frame_transport, name)
    for name in owner_names:
        assert getattr(electronic, name) is getattr(contract, name)
        assert getattr(frame_transport, name) is getattr(contract, name)
