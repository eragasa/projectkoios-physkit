"""Physical identities for closed-loop polar transport of electronic frames.

Optimization owns every unit-free transport calculation.  This module binds that
calculation to an exact ordered reduced-point route or a centered finite-periodic mesh,
an ordered ambient basis, a retained-frame identity, and an endpoint-sewing identity.
Frames, sewing maps, right-acting rotations, closure actions, endpoint mismatches,
correction steps, and closure defects are all dimensionless.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Self, final

import numpy as np
from numpy.typing import NDArray

from projectkoios.numerics.finite_periodic.meshes import CenteredUniformMesh
from projectkoios.optimization.alignment.closed_loop_polar_transport import (
    ClosedLoopPolarFrameTransportBackendFailure,
    ClosedLoopPolarFrameTransportBackendFailureReason,
    ClosedLoopPolarFrameTransportBackendOperation,
    ClosedLoopPolarFrameTransportLinkDiagnostic,
    ClosedLoopPolarFrameTransportLinkKind,
    ClosedLoopPolarFrameTransportLinkStatus,
    ClosedLoopPolarFrameTransportPrincipalPhaseBranch,
    ClosedLoopPolarFrameTransportStatus,
)
from projectkoios.optimization.alignment.numpy.closed_loop_polar_transport import (
    NumpyClosedLoopPolarFrameTransportDefinition,
    NumpyClosedLoopPolarFrameTransportResult,
)
from projectkoios.optimization.alignment.numpy.procrustes import (
    NumpyComplexUnitaryProcrustesSolver,
)
from projectkoios.optimization.alignment.procrustes import (
    DEFAULT_COMPLEX_UNITARY_PROCRUSTES_TOLERANCE,
    ComplexUnitaryProcrustesTolerance,
)
from projectkoios.optimization.alignment.scipy.closed_loop_polar_transport import (
    ScipyClosedLoopPolarFrameTransportSolver,
)
from projectkoios.physkit.core.results import ResultsObject

type ComplexArray = NDArray[np.complex128]
type ReducedPoint = tuple[float, ...]
type ReducedPointRoute = tuple[ReducedPoint, ...]


class FrameRotationAction(StrEnum):
    """Identify how retained-frame rotations act on frame matrices."""

    RIGHT = "right"


class EndpointSewingAction(StrEnum):
    """Identify how endpoint sewing acts in the ordered ambient basis."""

    LEFT_AMBIENT = "left-ambient"


def _identity(value: str, name: str) -> None:
    if type(value) is not str:
        raise TypeError(f"{name} must be a string")
    if not value:
        raise ValueError(f"{name} must be nonempty")


def _validate_identities(
    ambient_basis_identifier: str,
    retained_frame_identifier: str,
    endpoint_sewing_identifier: str,
) -> None:
    _identity(ambient_basis_identifier, "ambient_basis_identifier")
    _identity(retained_frame_identifier, "retained_frame_identifier")
    _identity(endpoint_sewing_identifier, "endpoint_sewing_identifier")


def _validate_reduced_points(
    reduced_points: ReducedPointRoute,
    point_count: int,
) -> None:
    if type(reduced_points) is not tuple:
        raise TypeError("reduced_points must be an exact tuple")
    if len(reduced_points) != point_count:
        raise ValueError("reduced_points length must equal the frame count")
    if not reduced_points:
        raise ValueError("reduced_points must be nonempty")
    dimension: int | None = None
    for point in reduced_points:
        if type(point) is not tuple:
            raise TypeError("each reduced point must be an exact tuple")
        if not point:
            raise ValueError("reduced points must have positive dimension")
        if dimension is None:
            dimension = len(point)
        elif len(point) != dimension:
            raise ValueError("all reduced points must have the same dimension")
        for component in point:
            if type(component) is not float:
                raise TypeError("reduced-point components must be built-in floats")
            if not math.isfinite(component):
                raise ValueError("reduced-point components must be finite")


def _numerical_definition(
    frames: ComplexArray,
    endpoint_sewing_map: ComplexArray,
    principal_phase_branch: ClosedLoopPolarFrameTransportPrincipalPhaseBranch,
    tolerance: ComplexUnitaryProcrustesTolerance,
) -> NumpyClosedLoopPolarFrameTransportDefinition:
    return NumpyClosedLoopPolarFrameTransportDefinition.from_immutable_arrays(
        frames,
        endpoint_sewing_map,
        principal_phase_branch,
        tolerance,
    )


@final
@dataclass(frozen=True, slots=True, eq=False)
class OrderedReducedPointFrameInventory:
    """Bind dimensionless frames to one exact ordered reduced-point route.

    Frame columns use ``ambient_basis_identifier`` and jointly represent the declared
    ``retained_frame_identifier``.  Local gauge rotations act on the right.  The
    endpoint sewing map acts on the left of the first frame in the ambient basis.
    The principal-phase branch is explicit. No coordinate support, tolerance,
    qualification, or acceptance policy is derived.
    """

    reduced_points: ReducedPointRoute
    frames: ComplexArray = field(repr=False)
    endpoint_sewing_map: ComplexArray = field(repr=False)
    ambient_basis_identifier: str
    retained_frame_identifier: str
    endpoint_sewing_identifier: str
    principal_phase_branch: ClosedLoopPolarFrameTransportPrincipalPhaseBranch
    tolerance: ComplexUnitaryProcrustesTolerance = (
        DEFAULT_COMPLEX_UNITARY_PROCRUSTES_TOLERANCE
    )

    @classmethod
    def from_arrays(
        cls,
        *,
        reduced_points: ReducedPointRoute,
        frames: ComplexArray,
        endpoint_sewing_map: ComplexArray,
        ambient_basis_identifier: str,
        retained_frame_identifier: str,
        endpoint_sewing_identifier: str,
        principal_phase_branch: ClosedLoopPolarFrameTransportPrincipalPhaseBranch,
        tolerance: ComplexUnitaryProcrustesTolerance = (
            DEFAULT_COMPLEX_UNITARY_PROCRUSTES_TOLERANCE
        ),
    ) -> Self:
        """Construct with immutable defensive copies owned by Optimization."""
        numerical = NumpyClosedLoopPolarFrameTransportDefinition.from_arrays(
            frames,
            endpoint_sewing_map,
            principal_phase_branch,
            tolerance,
        )
        return cls(
            reduced_points=reduced_points,
            frames=numerical.frames,
            endpoint_sewing_map=numerical.endpoint_sewing_map,
            ambient_basis_identifier=ambient_basis_identifier,
            retained_frame_identifier=retained_frame_identifier,
            endpoint_sewing_identifier=endpoint_sewing_identifier,
            principal_phase_branch=numerical.principal_phase_branch,
            tolerance=numerical.tolerance,
        )

    def __post_init__(self) -> None:
        numerical = self._numerical_definition()
        _validate_reduced_points(self.reduced_points, numerical.point_count)
        _validate_identities(
            self.ambient_basis_identifier,
            self.retained_frame_identifier,
            self.endpoint_sewing_identifier,
        )

    @property
    def frame_rotation_action(self) -> FrameRotationAction:
        """Return the fixed right-action convention ``frame @ rotation``."""
        return FrameRotationAction.RIGHT

    @property
    def endpoint_sewing_action(self) -> EndpointSewingAction:
        """Return the fixed convention ``sewing_map @ first_frame``."""
        return EndpointSewingAction.LEFT_AMBIENT

    @property
    def point_count(self) -> int:
        return int(self.frames.shape[0])

    @property
    def ambient_dimension(self) -> int:
        return int(self.frames.shape[1])

    @property
    def rank(self) -> int:
        return int(self.frames.shape[2])

    def _numerical_definition(self) -> NumpyClosedLoopPolarFrameTransportDefinition:
        return _numerical_definition(
            self.frames,
            self.endpoint_sewing_map,
            self.principal_phase_branch,
            self.tolerance,
        )


@final
@dataclass(frozen=True, slots=True, eq=False)
class FinitePeriodicMeshFrameInventory:
    """Bind dimensionless frames to one exact one-dimensional periodic mesh.

    The leading frame axis follows the one-dimensional mesh inventory in index order.
    The exact mesh object and explicit principal-phase branch are retained. Frame
    rotation and endpoint sewing use the same right and left-ambient conventions as
    :class:`OrderedReducedPointFrameInventory`.  Higher-dimensional meshes do not
    define one closed route or its independent cycle sewings and are rejected.
    """

    mesh: CenteredUniformMesh
    frames: ComplexArray = field(repr=False)
    endpoint_sewing_map: ComplexArray = field(repr=False)
    ambient_basis_identifier: str
    retained_frame_identifier: str
    endpoint_sewing_identifier: str
    principal_phase_branch: ClosedLoopPolarFrameTransportPrincipalPhaseBranch
    tolerance: ComplexUnitaryProcrustesTolerance = (
        DEFAULT_COMPLEX_UNITARY_PROCRUSTES_TOLERANCE
    )

    @classmethod
    def from_arrays(
        cls,
        *,
        mesh: CenteredUniformMesh,
        frames: ComplexArray,
        endpoint_sewing_map: ComplexArray,
        ambient_basis_identifier: str,
        retained_frame_identifier: str,
        endpoint_sewing_identifier: str,
        principal_phase_branch: ClosedLoopPolarFrameTransportPrincipalPhaseBranch,
        tolerance: ComplexUnitaryProcrustesTolerance = (
            DEFAULT_COMPLEX_UNITARY_PROCRUSTES_TOLERANCE
        ),
    ) -> Self:
        """Construct with immutable defensive copies owned by Optimization."""
        numerical = NumpyClosedLoopPolarFrameTransportDefinition.from_arrays(
            frames,
            endpoint_sewing_map,
            principal_phase_branch,
            tolerance,
        )
        return cls(
            mesh=mesh,
            frames=numerical.frames,
            endpoint_sewing_map=numerical.endpoint_sewing_map,
            ambient_basis_identifier=ambient_basis_identifier,
            retained_frame_identifier=retained_frame_identifier,
            endpoint_sewing_identifier=endpoint_sewing_identifier,
            principal_phase_branch=numerical.principal_phase_branch,
            tolerance=numerical.tolerance,
        )

    def __post_init__(self) -> None:
        if type(self.mesh) is not CenteredUniformMesh:
            raise TypeError("mesh must be CenteredUniformMesh")
        if self.mesh.dimension != 1:
            raise ValueError("mesh dimension must equal one for a closed route")
        numerical = self._numerical_definition()
        if numerical.point_count != self.mesh.point_count:
            raise ValueError(
                "frame count must equal the finite-periodic mesh point count"
            )
        _validate_identities(
            self.ambient_basis_identifier,
            self.retained_frame_identifier,
            self.endpoint_sewing_identifier,
        )

    @property
    def frame_rotation_action(self) -> FrameRotationAction:
        """Return the fixed right-action convention ``frame @ rotation``."""
        return FrameRotationAction.RIGHT

    @property
    def endpoint_sewing_action(self) -> EndpointSewingAction:
        """Return the fixed convention ``sewing_map @ first_frame``."""
        return EndpointSewingAction.LEFT_AMBIENT

    @property
    def point_count(self) -> int:
        return int(self.frames.shape[0])

    @property
    def ambient_dimension(self) -> int:
        return int(self.frames.shape[1])

    @property
    def rank(self) -> int:
        return int(self.frames.shape[2])

    def _numerical_definition(self) -> NumpyClosedLoopPolarFrameTransportDefinition:
        return _numerical_definition(
            self.frames,
            self.endpoint_sewing_map,
            self.principal_phase_branch,
            self.tolerance,
        )


type ClosedLoopFrameInventory = (
    OrderedReducedPointFrameInventory | FinitePeriodicMeshFrameInventory
)


@final
@dataclass(frozen=True, slots=True, eq=False, init=False)
class ClosedLoopFrameTransportOutcome(ResultsObject):
    """Retain exact physical identities with the complete Optimization result."""

    source: ClosedLoopFrameInventory
    numerical_result: NumpyClosedLoopPolarFrameTransportResult

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        raise TypeError(
            "ClosedLoopFrameTransportOutcome must be created by its owning action"
        )

    @classmethod
    def _from_action(
        cls,
        source: ClosedLoopFrameInventory,
        numerical_result: NumpyClosedLoopPolarFrameTransportResult,
    ) -> Self:
        result = object.__new__(cls)
        object.__setattr__(result, "source", source)
        object.__setattr__(result, "numerical_result", numerical_result)
        result._validate()
        return result

    def _validate(self) -> None:
        if type(self.source) not in (
            OrderedReducedPointFrameInventory,
            FinitePeriodicMeshFrameInventory,
        ):
            raise TypeError("source must be a closed-loop frame inventory")
        if type(self.numerical_result) is not NumpyClosedLoopPolarFrameTransportResult:
            raise TypeError(
                "numerical_result must be NumpyClosedLoopPolarFrameTransportResult"
            )
        definition = self.numerical_result.definition
        if (
            definition.frames is not self.source.frames
            or definition.endpoint_sewing_map is not self.source.endpoint_sewing_map
            or definition.principal_phase_branch
            is not self.source.principal_phase_branch
            or definition.tolerance is not self.source.tolerance
        ):
            raise ValueError(
                "numerical result definition is not correlated with its source"
            )
        # The nominal NumPy result validates transported frames against its retained
        # definition and frame_right_actions. Do not repeat that numerical operation.

    @property
    def numerical_definition(self) -> NumpyClosedLoopPolarFrameTransportDefinition:
        return self.numerical_result.definition

    @property
    def point_count(self) -> int:
        return self.numerical_result.point_count

    @property
    def ambient_dimension(self) -> int:
        return self.numerical_result.ambient_dimension

    @property
    def rank(self) -> int:
        return self.numerical_result.rank

    @property
    def principal_phase_branch(
        self,
    ) -> ClosedLoopPolarFrameTransportPrincipalPhaseBranch:
        return self.numerical_definition.principal_phase_branch

    @property
    def tolerance(self) -> ComplexUnitaryProcrustesTolerance:
        return self.numerical_definition.tolerance

    @property
    def status(self) -> ClosedLoopPolarFrameTransportStatus:
        return self.numerical_result.status

    @property
    def link_diagnostics(
        self,
    ) -> tuple[ClosedLoopPolarFrameTransportLinkDiagnostic, ...]:
        return self.numerical_result.link_diagnostics

    @property
    def deficient_link_indices(self) -> tuple[int, ...]:
        return self.numerical_result.deficient_link_indices

    @property
    def unattempted_link_indices(self) -> tuple[int, ...]:
        return self.numerical_result.unattempted_link_indices

    @property
    def transported_frames(self) -> ComplexArray | None:
        return self.numerical_result.transported_frames

    @property
    def frame_right_actions(self) -> ComplexArray | None:
        return self.numerical_result.frame_right_actions

    @property
    def closure_alignment_right_action(self) -> ComplexArray | None:
        return self.numerical_result.closure_alignment_right_action

    @property
    def closure_endpoint_mismatch(self) -> ComplexArray | None:
        return self.numerical_result.closure_endpoint_mismatch

    @property
    def distributed_closure_correction_step(self) -> ComplexArray | None:
        return self.numerical_result.distributed_closure_correction_step

    @property
    def closure_schur_vectors(self) -> ComplexArray | None:
        return self.numerical_result.closure_schur_vectors

    @property
    def maximum_absolute_entry_closure_alignment_residual(self) -> float | None:
        return self.numerical_result.maximum_absolute_entry_closure_alignment_residual

    @property
    def maximum_absolute_entry_closure_alignment_unitarity_defect(
        self,
    ) -> float | None:
        result = self.numerical_result
        return result.maximum_absolute_entry_closure_alignment_unitarity_defect

    @property
    def ordered_closure_endpoint_mismatch_phases(self) -> tuple[float, ...] | None:
        return self.numerical_result.ordered_closure_endpoint_mismatch_phases

    @property
    def svd_backend_identifier(self) -> str:
        return self.numerical_result.svd_backend_identifier

    @property
    def svd_backend_version(self) -> str:
        return self.numerical_result.svd_backend_version

    @property
    def schur_backend_identifier(self) -> str | None:
        return self.numerical_result.schur_backend_identifier

    @property
    def schur_backend_version(self) -> str | None:
        return self.numerical_result.schur_backend_version

    @property
    def frame_rotation_action(self) -> FrameRotationAction:
        return self.source.frame_rotation_action

    @property
    def endpoint_sewing_action(self) -> EndpointSewingAction:
        return self.source.endpoint_sewing_action

    @property
    def ambient_basis_identifier(self) -> str:
        return self.source.ambient_basis_identifier

    @property
    def retained_frame_identifier(self) -> str:
        return self.source.retained_frame_identifier

    @property
    def endpoint_sewing_identifier(self) -> str:
        return self.source.endpoint_sewing_identifier


@final
class ClosedLoopPolarFrameTransporter:
    """Delegate one physical frame inventory to the exact Optimization solver."""

    __slots__ = ()

    def transport(
        self,
        source: ClosedLoopFrameInventory,
        /,
    ) -> ClosedLoopFrameTransportOutcome:
        """Transport once without retry, fallback, or exception translation."""
        if type(source) not in (
            OrderedReducedPointFrameInventory,
            FinitePeriodicMeshFrameInventory,
        ):
            raise TypeError("source must be a closed-loop frame inventory")
        numerical = ScipyClosedLoopPolarFrameTransportSolver(
            NumpyComplexUnitaryProcrustesSolver()
        ).solve(source._numerical_definition())
        return ClosedLoopFrameTransportOutcome._from_action(source, numerical)


__all__ = [
    "ClosedLoopFrameTransportOutcome",
    "ClosedLoopPolarFrameTransportBackendFailure",
    "ClosedLoopPolarFrameTransportBackendFailureReason",
    "ClosedLoopPolarFrameTransportBackendOperation",
    "ClosedLoopPolarFrameTransportLinkDiagnostic",
    "ClosedLoopPolarFrameTransportLinkKind",
    "ClosedLoopPolarFrameTransportLinkStatus",
    "ClosedLoopPolarFrameTransportPrincipalPhaseBranch",
    "ClosedLoopPolarFrameTransportStatus",
    "ClosedLoopPolarFrameTransporter",
    "EndpointSewingAction",
    "FinitePeriodicMeshFrameInventory",
    "FrameRotationAction",
    "OrderedReducedPointFrameInventory",
]
