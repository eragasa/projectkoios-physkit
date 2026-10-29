# Module `projectkoios.physkit.solidstate.electronic.frame_transport`

## Responsibility

This module is the physical-semantic boundary around Optimization's unit-free
closed-loop polar frame transport. It does not implement a second transport algorithm.
It binds each calculation to:

- either an exact ordered tuple of finite reduced-coordinate points or an exact
  one-dimensional `CenteredUniformMesh` object in index order;
- a nonempty ordered ambient-basis identity;
- a nonempty retained-frame identity;
- a nonempty endpoint-sewing identity;
- an explicit principal-phase branch, either `(-pi, pi]` or `[-pi, pi)`;
- right-acting retained-frame rotations, `frame @ rotation`; and
- left ambient endpoint sewing, `endpoint_sewing_map @ first_frame`.

Frames, sewing maps, rotations, closure actions, endpoint mismatches, correction
steps, singular spectra, cutoffs, margins, maximum-absolute-entry defects, and phases
are dimensionless. The wrapper does not attach an energy or length unit to them.

## Public records and action

- `OrderedReducedPointFrameInventory` owns immutable complex128 frames and endpoint
  sewing while preserving the supplied reduced-point tuple by identity.
- `FinitePeriodicMeshFrameInventory` owns the same arrays and preserves a supplied
  one-dimensional `CenteredUniformMesh` by identity. Its frame count must equal
  `mesh.point_count`. Two- and three-dimensional meshes are rejected because one
  flattened traversal does not define their independent periodic-cycle sewings.
- `FrameRotationAction.RIGHT` and `EndpointSewingAction.LEFT_AMBIENT` make the two
  matrix-action conventions explicit.
- `ClosedLoopPolarFrameTransporter.transport(source)` constructs the required
  `NumpyComplexUnitaryProcrustesSolver` collaborator and invokes exactly one
  `ScipyClosedLoopPolarFrameTransportSolver.solve(...)` call.
- `ClosedLoopFrameTransportOutcome` retains the exact source object and complete
  `NumpyClosedLoopPolarFrameTransportResult`. It correlates the source to the result's
  retained definition arrays by identity. Convenience properties expose status, link
  diagnostics, deficient and unattempted indices, transported frames, frame right
  actions, closure alignment right action, closure endpoint mismatch, distributed
  closure correction step, Schur vectors, storage-ordered mismatch phases, explicitly
  named maximum-absolute-entry defects, backend evidence, and physical identities.

The module and `projectkoios.physkit.solidstate.electronic` façade also export the exact
Optimization status, link-kind, link-status, link-diagnostic, principal-phase-branch,
backend-operation, backend-failure-reason, and backend-failure type identities. Backend
failures pass through unchanged with their operation, reason, and exact link index.

## Numerical boundary

Optimization's ownership is split across three exact layers:

- `projectkoios.optimization.alignment.closed_loop_polar_transport` owns dependency-free
  statuses, link diagnostics, and typed backend failures;
- `projectkoios.optimization.alignment.numpy.closed_loop_polar_transport` owns
  `NumpyClosedLoopPolarFrameTransportDefinition` and the self-validating
  `NumpyClosedLoopPolarFrameTransportResult`; and
- `projectkoios.optimization.alignment.scipy.closed_loop_polar_transport` owns only
  `ScipyClosedLoopPolarFrameTransportSolver` execution.

The PhysKit inventories use the NumPy definition's immutable complex128 storage, its
required `ClosedLoopPolarFrameTransportPrincipalPhaseBranch`, and its existing
`ComplexUnitaryProcrustesTolerance`; PhysKit introduces no numerical threshold. The
NumPy result retains the exact definition, immutable frame right actions, and the
correlation `transported_frames == definition.frames @ frame_right_actions` within its
owner validation bound.

For the closure protocol, Optimization reports:

- `closure_alignment_right_action` A, aligning the transported last frame to the sewn
  first frame;
- `closure_endpoint_mismatch` E = A† in KSD orientation;
- `distributed_closure_correction_step` W with W^M = E† = A; and
- `ordered_closure_endpoint_mismatch_phases`, whose ascending order is deterministic
  storage order only and carries no band labeling claim.

Every reported defect explicitly uses the maximum-absolute-entry metric. Optimization
remains responsible for sequential right-unitary Procrustes transport, rank
qualification, closure Schur decomposition, distributed closure, diagnostics, and
backend evidence.

Rank-deficient results are ordinary typed outcomes. The wrapper preserves the first
rank-deficient link and every explicitly unattempted dependent link. It does not turn
that status into an exception or manufacture transported frames. Typed SVD and Schur
backend failures remain exceptions and are not inferred by parsing text.

## Exclusions and evidence limits

This module chooses no point route, mesh, retained subspace, sewing map, tolerance,
support, KSD status, scientific threshold, protocol acceptance, or claim. Passing its
focused tests establishes only representation correlation, delegation, and the tested
numerical examples; it does not establish physical adequacy of a route or retained
frame.

## KSD consumption route

A KSD consumer constructs either inventory from its already selected dimensionless
frame matrices and sewing map, supplies its exact route or mesh and identity strings,
and explicitly selects
`ClosedLoopPolarFrameTransportPrincipalPhaseBranch.negative_pi_open` for the KSD
`(-pi, pi]` protocol before calling `ClosedLoopPolarFrameTransporter.transport`. It
branches on `ClosedLoopPolarFrameTransportStatus`, consumes typed link diagnostics,
backend evidence, immutable frame right actions, A, E, W, and mismatch phases from the
outcome, and retains `outcome.source` plus `outcome.numerical_definition` for exact
identity correlation. Any normalized Frobenius diagnostics remain KSD-side; PhysKit
only exposes Optimization's explicitly named maximum-absolute-entry evidence. KSD
remains responsible for all material-specific interpretation and acceptance.

## Local mapping

- Code: `src/python/projectkoios/physkit/solidstate/electronic/frame_transport.py`
- Tests: `tests/projectkoios/physkit/solidstate/electronic/frame_transport/`
