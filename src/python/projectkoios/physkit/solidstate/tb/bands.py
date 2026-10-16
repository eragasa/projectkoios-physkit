r"""Complete unit-aware bands from matrix-valued Bloch Hamiltonians.

The public adapters evaluate one :class:`BlockHoppingModel` or an explicit ordered
source-and-truncation family, apply a caller-supplied numerical Hermiticity tolerance,
construct the explicit Hermitian representative :math:`(H + H^\dagger)/2`, and return
every ascending eigenvalue at every supplied fractional reciprocal coordinate.

Shared-grid finite-Fourier evaluation and the generic unit-free spectrum operation are
delegated to pinned ``projectkoios.numerics`` owners. The family path shares phases,
concatenates matrices, and partitions results without selecting a truncation policy.
No scientific tolerance, gap metric, cutoff policy, acceptance decision, cyclic
displacement interpretation, or finite-mesh Nyquist rule is owned here.

Notes
-----
For each fractional reciprocal coordinate :math:`\mathbf q_p`, the Bloch adapter
provides a square matrix :math:`A_p` in the model's ordered basis. Accepted inputs
produce the Hermitian representative

.. math::

   H_p = \frac{A_p + A_p^\dagger}{2}

and one complete ascending spectrum :math:`E_{p,0}\leq\cdots\leq E_{p,n-1}`.
The returned energy matrix therefore has shape ``(point_count, basis_size)``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from projectkoios.numerics.linalg.hermitian_spectrum import (
    BatchedHermitianSpectrumDefinition,
    BatchedHermitianSpectrumTolerance,
)
from projectkoios.numerics.linalg.numpy.hermitian_spectrum import (
    NumpyBatchedHermitianSpectrumSolver,
)
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.lattice.finite_domain import LatticeDisplacement
from projectkoios.physkit.periodic.lattice.hopping import BlockHoppingModel
from projectkoios.physkit.solidstate.tb.bloch import (
    BlockHoppingBlochEvaluation,
    BlockHoppingBlochEvaluator,
)
from projectkoios.physkit.solidstate.tb.truncation import (
    BlockHoppingTruncationResult,
    ConjugateClosedBlockHoppingTruncator,
)
from projectkoios.physkit.units import MatrixQuantity, ScalarQuantity


@dataclass(frozen=True, slots=True)
class BlockHoppingBandTolerance:
    r"""Caller-supplied numerical Hermiticity tolerances.

    Parameters
    ----------
    hermiticity_absolute : ScalarQuantity
        Finite nonnegative absolute tolerance. Evaluation requires the hopping
        model's energy unit.
    hermiticity_relative : float
        Exact finite binary64 relative tolerance in the interval ``[0, 1)``.

    Notes
    -----
    For an input matrix :math:`A`, absolute tolerance :math:`a`, and relative
    tolerance :math:`r`, the numerical cutoff is

    .. math::

       c(A) = a + r\max_{ij}|A_{ij}|.

    This record carries numerical parameters only. It does not choose a scientific
    tolerance or establish model acceptance.
    """

    hermiticity_absolute: ScalarQuantity
    hermiticity_relative: float

    def __post_init__(self) -> None:
        """The initializer validates absolute and relative tolerance bounds."""
        if type(self.hermiticity_absolute) is not ScalarQuantity:
            raise TypeError("hermiticity_absolute must be ScalarQuantity")
        if self.hermiticity_absolute.magnitude < 0.0:
            raise ValueError("hermiticity_absolute must be nonnegative")
        if type(self.hermiticity_relative) is not float:
            raise TypeError("hermiticity_relative must be a built-in float")
        if not math.isfinite(self.hermiticity_relative):
            raise ValueError("hermiticity_relative must be finite")
        if not 0.0 <= self.hermiticity_relative < 1.0:
            raise ValueError("hermiticity_relative must lie in [0, 1)")


@dataclass(frozen=True, slots=True)
class BlockHoppingRangeBandSelection:
    """One explicit ordinary-lattice inventory in an ordered range family.

    Parameters
    ----------
    retained_displacements : tuple of LatticeDisplacement
        Nonempty, sorted, unique inventory closed under ordinary negation.
    truncated_identifier : str
        Nonempty identity for the retained model. Family evaluation additionally
        requires identifiers and inventories to be unique across selections.

    Notes
    -----
    The record contains no distance metric, cutoff inference, cyclic equivalence,
    Nyquist rule, error tolerance, or acceptance policy.
    """

    retained_displacements: tuple[LatticeDisplacement, ...]
    truncated_identifier: str

    def __post_init__(self) -> None:
        """Validate one explicit conjugate-closed ordinary-lattice selection."""
        if type(self.retained_displacements) is not tuple or any(
            type(displacement) is not LatticeDisplacement
            for displacement in self.retained_displacements
        ):
            raise TypeError(
                "retained_displacements must be a tuple of LatticeDisplacement values"
            )
        if not self.retained_displacements:
            raise ValueError("retained_displacements must be nonempty")
        dimension = self.retained_displacements[0].dimension
        if any(
            displacement.dimension is not dimension
            for displacement in self.retained_displacements[1:]
        ):
            raise ValueError("retained displacement dimensions must match")
        components = tuple(
            displacement.components for displacement in self.retained_displacements
        )
        if components != tuple(sorted(set(components))):
            raise ValueError("retained_displacements must be sorted and unique")
        inventory = set(components)
        if any(
            tuple(-component for component in displacement) not in inventory
            for displacement in components
        ):
            raise ValueError("retained_displacements must be closed under negation")
        if type(self.truncated_identifier) is not str:
            raise TypeError("truncated_identifier must be a string")
        if not self.truncated_identifier:
            raise ValueError("truncated_identifier must be nonempty")


@dataclass(frozen=True, slots=True)
class BlockHoppingBandPointDiagnostics:
    r"""Raw numerical facts for one Bloch-matrix input row.

    Parameters
    ----------
    point_index : int
        Zero-based row index in the retained fractional-coordinate matrix.
    input_hermiticity_maximum_absolute_defect : ScalarQuantity
        Maximum entrywise magnitude of :math:`A-A^\dagger`.
    representative_adjustment_maximum_absolute_defect : ScalarQuantity
        Maximum entrywise magnitude of :math:`(A+A^\dagger)/2-A`.
    input_matrix_maximum_absolute_entry : ScalarQuantity
        Maximum entrywise magnitude of the input matrix.
    applied_hermiticity_cutoff : ScalarQuantity
        Numerical cutoff applied to the input Hermiticity defect.

    Notes
    -----
    Every diagnostic quantity carries the model energy unit. A successful point
    satisfies ``input_hermiticity_maximum_absolute_defect <=
    applied_hermiticity_cutoff``; equality is accepted.
    """

    point_index: int
    input_hermiticity_maximum_absolute_defect: ScalarQuantity
    representative_adjustment_maximum_absolute_defect: ScalarQuantity
    input_matrix_maximum_absolute_entry: ScalarQuantity
    applied_hermiticity_cutoff: ScalarQuantity

    def __post_init__(self) -> None:
        """The initializer validates index, unit, magnitude, and cutoff invariants."""
        if type(self.point_index) is not int:
            raise TypeError("point_index must be an exact built-in integer")
        if self.point_index < 0:
            raise ValueError("point_index must be nonnegative")
        quantities = (
            self.input_hermiticity_maximum_absolute_defect,
            self.representative_adjustment_maximum_absolute_defect,
            self.input_matrix_maximum_absolute_entry,
            self.applied_hermiticity_cutoff,
        )
        if any(type(quantity) is not ScalarQuantity for quantity in quantities):
            raise TypeError("point diagnostic values must be ScalarQuantity values")
        if any(quantity.magnitude < 0.0 for quantity in quantities):
            raise ValueError("point diagnostic values must be nonnegative")
        if any(quantity.unit != quantities[0].unit for quantity in quantities[1:]):
            raise ValueError("point diagnostic values must use one common unit")
        if (
            self.input_hermiticity_maximum_absolute_defect.magnitude
            > self.applied_hermiticity_cutoff.magnitude
        ):
            raise ValueError("Hermiticity defect must not exceed the applied cutoff")


@dataclass(frozen=True, slots=True)
class BlockHoppingBandDiagnostics:
    """Componentwise maxima across a nonempty Bloch-matrix batch.

    Parameters
    ----------
    maximum_input_hermiticity_absolute_defect : ScalarQuantity
        Largest per-point input Hermiticity defect.
    maximum_representative_adjustment_absolute_defect : ScalarQuantity
        Largest per-point representative adjustment.
    maximum_input_matrix_absolute_entry : ScalarQuantity
        Largest per-point input maximum-absolute entry.
    maximum_applied_hermiticity_cutoff : ScalarQuantity
        Largest per-point applied Hermiticity cutoff.

    Notes
    -----
    Empty band evaluations use ``diagnostics=None`` because maxima over an empty
    batch are undefined.
    """

    maximum_input_hermiticity_absolute_defect: ScalarQuantity
    maximum_representative_adjustment_absolute_defect: ScalarQuantity
    maximum_input_matrix_absolute_entry: ScalarQuantity
    maximum_applied_hermiticity_cutoff: ScalarQuantity

    def __post_init__(self) -> None:
        """The initializer validates aggregate magnitude and unit invariants."""
        quantities = (
            self.maximum_input_hermiticity_absolute_defect,
            self.maximum_representative_adjustment_absolute_defect,
            self.maximum_input_matrix_absolute_entry,
            self.maximum_applied_hermiticity_cutoff,
        )
        if any(type(quantity) is not ScalarQuantity for quantity in quantities):
            raise TypeError("aggregate diagnostic values must be ScalarQuantity values")
        if any(quantity.magnitude < 0.0 for quantity in quantities):
            raise ValueError("aggregate diagnostic values must be nonnegative")
        if any(quantity.unit != quantities[0].unit for quantity in quantities[1:]):
            raise ValueError("aggregate diagnostic values must use one common unit")
        if (
            self.maximum_input_hermiticity_absolute_defect.magnitude
            > self.maximum_applied_hermiticity_cutoff.magnitude
        ):
            raise ValueError(
                "maximum Hermiticity defect must not exceed the maximum cutoff"
            )


def _aggregate_diagnostics(
    point_diagnostics: tuple[BlockHoppingBandPointDiagnostics, ...],
    model: BlockHoppingModel,
) -> BlockHoppingBandDiagnostics | None:
    """The reduction returns componentwise maxima for one band evaluation.

    Parameters
    ----------
    point_diagnostics : tuple of BlockHoppingBandPointDiagnostics
        Ordered point diagnostics to reduce.
    model : BlockHoppingModel
        Model supplying the energy unit.

    Returns
    -------
    BlockHoppingBandDiagnostics or None
        Componentwise maxima, or ``None`` when no points exist.
    """
    if not point_diagnostics:
        return None

    def maximum(attribute: str) -> ScalarQuantity:
        return ScalarQuantity(
            max(
                (
                    getattr(diagnostic, attribute).magnitude
                    for diagnostic in point_diagnostics
                ),
            ),
            model.energy_unit,
        )

    return BlockHoppingBandDiagnostics(
        maximum_input_hermiticity_absolute_defect=maximum(
            "input_hermiticity_maximum_absolute_defect"
        ),
        maximum_representative_adjustment_absolute_defect=maximum(
            "representative_adjustment_maximum_absolute_defect"
        ),
        maximum_input_matrix_absolute_entry=maximum(
            "input_matrix_maximum_absolute_entry"
        ),
        maximum_applied_hermiticity_cutoff=maximum("applied_hermiticity_cutoff"),
    )


@dataclass(frozen=True, slots=True, eq=False)
class BlockHoppingBandEvaluation(ResultsObject):
    """A correlated Bloch evaluation, complete bands, and diagnostics.

    Parameters
    ----------
    bloch_evaluation : BlockHoppingBlochEvaluation
        Authoritative Bloch matrices, model, and fractional coordinates used by the
        spectrum calculation.
    band_energies : MatrixQuantity
        Complete ascending spectra with shape ``(point_count, basis_size)`` in the
        model energy unit.
    tolerance : BlockHoppingBandTolerance
        Caller-supplied tolerance used for the whole batch.
    point_diagnostics : tuple of BlockHoppingBandPointDiagnostics
        Per-point diagnostics in fractional-coordinate row order.
    diagnostics : BlockHoppingBandDiagnostics or None
        Componentwise maxima for a nonempty batch, or ``None`` for an empty batch.
    backend_identifier : str
        Nonempty identifier for the selected complete-spectrum backend.
    backend_version : str
        Nonempty backend version evidence.

    Notes
    -----
    Row ``p`` of ``band_energies`` and element ``p`` of ``point_diagnostics`` both
    refer to row ``p`` of ``bloch_evaluation.fractional_wave_vectors`` and matrix
    ``p`` of ``bloch_evaluation.bloch_matrices``. Numerical derivation is owned by
    :class:`BlockHoppingBandEvaluator`, avoiding a second eigensolver call.
    """

    bloch_evaluation: BlockHoppingBlochEvaluation
    band_energies: MatrixQuantity
    tolerance: BlockHoppingBandTolerance
    point_diagnostics: tuple[BlockHoppingBandPointDiagnostics, ...]
    diagnostics: BlockHoppingBandDiagnostics | None
    backend_identifier: str
    backend_version: str

    def __post_init__(self) -> None:
        """The initializer validates structural correlation and evidence shapes."""
        if type(self.bloch_evaluation) is not BlockHoppingBlochEvaluation:
            raise TypeError("bloch_evaluation must be BlockHoppingBlochEvaluation")
        point_count = len(self.bloch_evaluation.bloch_matrices)
        # One row of energies corresponds to one retained Bloch matrix and q point.
        if type(self.band_energies) is not MatrixQuantity:
            raise TypeError("band_energies must be MatrixQuantity")
        if self.band_energies.unit != self.model.energy_unit:
            raise ValueError("band_energies must use the model energy unit")
        if self.band_energies.magnitude.shape != (
            point_count,
            self.model.basis_size,
        ):
            raise ValueError("band_energies shape must be (point_count, basis_size)")
        # Adjacent-column comparison verifies ascending order without re-solving.
        if np.any(
            self.band_energies.magnitude[:, 1:] < self.band_energies.magnitude[:, :-1]
        ):
            raise ValueError("band_energies rows must be ascending")
        if type(self.tolerance) is not BlockHoppingBandTolerance:
            raise TypeError("tolerance must be BlockHoppingBandTolerance")
        if self.tolerance.hermiticity_absolute.unit != self.model.energy_unit:
            raise ValueError("absolute tolerance must use the model energy unit")
        if type(self.point_diagnostics) is not tuple or any(
            type(diagnostic) is not BlockHoppingBandPointDiagnostics
            for diagnostic in self.point_diagnostics
        ):
            raise TypeError(
                "point_diagnostics must contain BlockHoppingBandPointDiagnostics"
            )
        if len(self.point_diagnostics) != point_count:
            raise ValueError("point_diagnostics must match the Bloch-matrix count")
        for index, diagnostic in enumerate(self.point_diagnostics):
            if diagnostic.point_index != index:
                raise ValueError("point_diagnostics must use canonical point indices")
            if (
                diagnostic.input_hermiticity_maximum_absolute_defect.unit
                != self.model.energy_unit
            ):
                raise ValueError("point diagnostics must use the model energy unit")
            # The reported cutoff remains reproducible from retained public inputs.
            expected_cutoff = (
                self.tolerance.hermiticity_absolute.magnitude
                + self.tolerance.hermiticity_relative
                * diagnostic.input_matrix_maximum_absolute_entry.magnitude
            )
            if diagnostic.applied_hermiticity_cutoff.magnitude != expected_cutoff:
                raise ValueError(
                    "applied cutoffs must match the supplied tolerance and input scale"
                )
        if (
            self.diagnostics is not None
            and type(self.diagnostics) is not BlockHoppingBandDiagnostics
        ):
            raise TypeError("diagnostics must be BlockHoppingBandDiagnostics or None")
        if self.diagnostics != _aggregate_diagnostics(
            self.point_diagnostics,
            self.model,
        ):
            raise ValueError("diagnostics must equal componentwise point maxima")
        if type(self.backend_identifier) is not str:
            raise TypeError("backend_identifier must be a string")
        if not self.backend_identifier:
            raise ValueError("backend_identifier must be nonempty")
        if type(self.backend_version) is not str:
            raise TypeError("backend_version must be a string")
        if not self.backend_version:
            raise ValueError("backend_version must be nonempty")

    @property
    def model(self) -> BlockHoppingModel:
        """The property returns the authoritative Bloch evaluation's model."""
        return self.bloch_evaluation.model

    @property
    def fractional_wave_vectors(self) -> MatrixQuantity:
        """The property returns the authoritative fractional coordinates."""
        return self.bloch_evaluation.fractional_wave_vectors


@dataclass(frozen=True, slots=True, eq=False)
class BlockHoppingRangeBandEntry(ResultsObject):
    """Bind one exact truncation result to its exact band evaluation.

    Parameters
    ----------
    truncation : BlockHoppingTruncationResult
        Metadata-preserving explicit source-term selection and omitted power.
    band_evaluation : BlockHoppingBandEvaluation
        Complete bands of the exact retained model.
    """

    truncation: BlockHoppingTruncationResult
    band_evaluation: BlockHoppingBandEvaluation

    def __post_init__(self) -> None:
        """Validate exact retained-model identity across both results."""
        if type(self.truncation) is not BlockHoppingTruncationResult:
            raise TypeError("truncation must be BlockHoppingTruncationResult")
        if type(self.band_evaluation) is not BlockHoppingBandEvaluation:
            raise TypeError("band_evaluation must be BlockHoppingBandEvaluation")
        if self.band_evaluation.model is not self.truncation.truncated:
            raise ValueError("band evaluation must retain the exact truncated model")


@dataclass(frozen=True, slots=True, eq=False)
class BlockHoppingRangeBandEvaluation(ResultsObject):
    """Correlate source bands with an ordered explicit truncation family.

    Parameters
    ----------
    source_band_evaluation : BlockHoppingBandEvaluation
        Complete source-model matrices, bands, diagnostics, and backend evidence.
    entries : tuple of BlockHoppingRangeBandEntry
        Nonempty ordered truncation and band pairs.

    Notes
    -----
    All evaluations retain the same coordinate and tolerance objects. Inventories and
    retained-model identifiers are unique, but no numerical range or cutoff policy is
    inferred from their order.
    """

    source_band_evaluation: BlockHoppingBandEvaluation
    entries: tuple[BlockHoppingRangeBandEntry, ...]

    def __post_init__(self) -> None:
        """Validate source, ordering, uniqueness, and shared evaluation inputs."""
        if type(self.source_band_evaluation) is not BlockHoppingBandEvaluation:
            raise TypeError("source_band_evaluation must be BlockHoppingBandEvaluation")
        if type(self.entries) is not tuple or any(
            type(entry) is not BlockHoppingRangeBandEntry for entry in self.entries
        ):
            raise TypeError("entries must be a tuple of BlockHoppingRangeBandEntry")
        if not self.entries:
            raise ValueError("entries must be nonempty")
        source = self.source_band_evaluation.model
        coordinates = self.source_band_evaluation.fractional_wave_vectors
        tolerance = self.source_band_evaluation.tolerance
        for entry in self.entries:
            if entry.truncation.source is not source:
                raise ValueError("every truncation must retain the exact source model")
            if entry.band_evaluation.fractional_wave_vectors is not coordinates:
                raise ValueError("every band evaluation must retain shared coordinates")
            if entry.band_evaluation.tolerance is not tolerance:
                raise ValueError(
                    "every band evaluation must retain the shared tolerance"
                )
        identifiers = tuple(
            entry.truncation.truncated.identifier for entry in self.entries
        )
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("truncated model identifiers must be unique")
        inventories = tuple(
            tuple(
                displacement.components
                for displacement in entry.truncation.retained_displacements
            )
            for entry in self.entries
        )
        if len(set(inventories)) != len(inventories):
            raise ValueError("retained displacement inventories must be unique")

    @property
    def source(self) -> BlockHoppingModel:
        """Return the exact source model retained by the source band result."""
        return self.source_band_evaluation.model

    @property
    def fractional_wave_vectors(self) -> MatrixQuantity:
        """Return the exact shared reciprocal-coordinate record."""
        return self.source_band_evaluation.fractional_wave_vectors

    @property
    def tolerance(self) -> BlockHoppingBandTolerance:
        """Return the exact shared Hermiticity tolerance."""
        return self.source_band_evaluation.tolerance


class BlockHoppingBandEvaluator:
    """Complete block-hopping bands under caller-supplied tolerances.

    Parameters
    ----------
    bloch_evaluator : BlockHoppingBlochEvaluator or None, optional
        Matrix Bloch evaluator to consume. ``None`` constructs the maintained
        default evaluator.

    Notes
    -----
    The evaluator owns only the unit-aware physics adaptation. The pinned Numerics
    owner performs generic Hermiticity validation, representative construction, and
    complete-spectrum calculation.

    Each nonempty evaluation follows four numerical stages: Bloch-matrix assembly,
    maximum-absolute Hermiticity diagnostics, explicit Hermitian representative
    construction, and one batched complete-spectrum solve. Empty evaluations retain
    their matrix dimension without invoking linear algebra.
    """

    __slots__ = ("_bloch_evaluator", "_spectrum_solver")

    def __init__(
        self,
        bloch_evaluator: BlockHoppingBlochEvaluator | None = None,
    ) -> None:
        """The constructor retains a Bloch evaluator and canonical backend."""
        if bloch_evaluator is None:
            bloch_evaluator = BlockHoppingBlochEvaluator()
        elif not isinstance(bloch_evaluator, BlockHoppingBlochEvaluator):
            raise TypeError(
                "bloch_evaluator must be BlockHoppingBlochEvaluator or None"
            )
        self._bloch_evaluator = bloch_evaluator
        self._spectrum_solver = NumpyBatchedHermitianSpectrumSolver()

    @property
    def bloch_evaluator(self) -> BlockHoppingBlochEvaluator:
        """The property returns the matrix Bloch evaluator used by this adapter.

        Returns
        -------
        BlockHoppingBlochEvaluator
            Exact nominal evaluator retained at construction.
        """
        return self._bloch_evaluator

    def evaluate_bloch(
        self,
        bloch_evaluation: BlockHoppingBlochEvaluation,
        tolerance: BlockHoppingBandTolerance,
    ) -> BlockHoppingBandEvaluation:
        """Evaluate complete bands from one existing Bloch evaluation.

        Parameters
        ----------
        bloch_evaluation : BlockHoppingBlochEvaluation
            Exact model, coordinates, and Bloch matrices to retain in the result.
        tolerance : BlockHoppingBandTolerance
            Caller-supplied numerical Hermiticity tolerance.

        Returns
        -------
        BlockHoppingBandEvaluation
            Complete ascending bands and correlated numerical diagnostics.

        Raises
        ------
        TypeError
            If an argument has unsupported nominal ownership.
        ValueError
            If units are incompatible or a Hermiticity defect exceeds its cutoff.
        FloatingPointError
            If diagnostic, representative, or spectral arithmetic becomes nonfinite.
        RuntimeError
            If the complete Hermitian eigensolver does not converge.

        Notes
        -----
        The supplied Bloch evaluation is retained exactly and is not recomputed.
        """
        if type(bloch_evaluation) is not BlockHoppingBlochEvaluation:
            raise TypeError("bloch_evaluation must be BlockHoppingBlochEvaluation")
        return self.evaluate_bloch_family((bloch_evaluation,), tolerance)[0]

    def evaluate_bloch_family(
        self,
        bloch_evaluations: tuple[BlockHoppingBlochEvaluation, ...],
        tolerance: BlockHoppingBandTolerance,
    ) -> tuple[BlockHoppingBandEvaluation, ...]:
        """Evaluate a compatible Bloch family with one complete-spectrum batch.

        Parameters
        ----------
        bloch_evaluations : tuple of BlockHoppingBlochEvaluation
            Nonempty ordered evaluations sharing coordinates, basis metadata, energy
            semantics, and point count.
        tolerance : BlockHoppingBandTolerance
            Numerical Hermiticity tolerance applied independently to every matrix.

        Returns
        -------
        tuple of BlockHoppingBandEvaluation
            Results partitioned in supplied evaluation order.

        Raises
        ------
        TypeError
            If inputs do not have their exact nominal PhysKit types.
        ValueError
            If family metadata, coordinates, units, or Hermiticity gates disagree.
        FloatingPointError
            If diagnostic, representative, or eigensolver output arithmetic is
            nonfinite.
        RuntimeError
            If the complete Hermitian eigensolver does not converge.

        Notes
        -----
        All nonempty matrices are concatenated into one canonical Numerics solve and
        partitioned deterministically. Empty shared coordinates skip linear algebra.
        """
        if type(bloch_evaluations) is not tuple or any(
            type(evaluation) is not BlockHoppingBlochEvaluation
            for evaluation in bloch_evaluations
        ):
            raise TypeError(
                "bloch_evaluations must be a tuple of BlockHoppingBlochEvaluation"
            )
        if not bloch_evaluations:
            raise ValueError("bloch_evaluations must be nonempty")
        if type(tolerance) is not BlockHoppingBandTolerance:
            raise TypeError("tolerance must be BlockHoppingBandTolerance")
        reference = bloch_evaluations[0]
        model = reference.model
        coordinates = reference.fractional_wave_vectors
        point_count = len(reference.bloch_matrices)
        if tolerance.hermiticity_absolute.unit != model.energy_unit:
            raise ValueError("absolute tolerance must use the model energy unit")
        for evaluation in bloch_evaluations[1:]:
            family_model = evaluation.model
            if evaluation.fractional_wave_vectors is not coordinates:
                raise ValueError("Bloch family must retain exact shared coordinates")
            if family_model.dimension is not model.dimension:
                raise ValueError("Bloch family model dimensions must match")
            if family_model.basis_size != model.basis_size:
                raise ValueError("Bloch family basis sizes must match")
            if family_model.basis_identifier != model.basis_identifier:
                raise ValueError("Bloch family basis identifiers must match")
            if family_model.energy_unit != model.energy_unit:
                raise ValueError("Bloch family energy units must match")
            if family_model.energy_reference != model.energy_reference:
                raise ValueError("Bloch family energy references must match")
            if len(evaluation.bloch_matrices) != point_count:
                raise ValueError("Bloch family point counts must match")

        # Evaluation-major concatenation makes each model occupy one contiguous block
        # of point-major matrices for deterministic result partitioning.
        matrices = tuple(
            tuple(tuple(complex(entry) for entry in row) for row in matrix.magnitude)
            for evaluation in bloch_evaluations
            for matrix in evaluation.bloch_matrices
        )
        numerical_tolerance = BatchedHermitianSpectrumTolerance(
            tolerance.hermiticity_absolute.magnitude,
            tolerance.hermiticity_relative,
        )
        if matrices:
            definition = BatchedHermitianSpectrumDefinition(
                matrices,
                numerical_tolerance,
            )
        else:
            definition = BatchedHermitianSpectrumDefinition.empty(
                model.basis_size,
                numerical_tolerance,
            )
        spectrum = self._spectrum_solver.solve(definition)

        results: list[BlockHoppingBandEvaluation] = []
        for family_index, evaluation in enumerate(bloch_evaluations):
            start = family_index * point_count
            stop = start + point_count
            model = evaluation.model
            owner_diagnostics = spectrum.matrix_diagnostics[start:stop]
            point_diagnostics = tuple(
                BlockHoppingBandPointDiagnostics(
                    point_index=index,
                    input_hermiticity_maximum_absolute_defect=ScalarQuantity(
                        diagnostic.input_hermiticity_maximum_absolute_defect,
                        model.energy_unit,
                    ),
                    representative_adjustment_maximum_absolute_defect=ScalarQuantity(
                        diagnostic.representative_adjustment_maximum_absolute_defect,
                        model.energy_unit,
                    ),
                    input_matrix_maximum_absolute_entry=ScalarQuantity(
                        diagnostic.input_matrix_maximum_absolute_entry,
                        model.energy_unit,
                    ),
                    applied_hermiticity_cutoff=ScalarQuantity(
                        diagnostic.applied_hermiticity_cutoff,
                        model.energy_unit,
                    ),
                )
                for index, diagnostic in enumerate(owner_diagnostics)
            )
            eigenvalues = spectrum.eigenvalues[start:stop]
            results.append(
                BlockHoppingBandEvaluation(
                    bloch_evaluation=evaluation,
                    band_energies=MatrixQuantity(
                        np.asarray(eigenvalues, dtype=np.float64).reshape(
                            point_count,
                            model.basis_size,
                        ),
                        model.energy_unit,
                    ),
                    tolerance=tolerance,
                    point_diagnostics=point_diagnostics,
                    diagnostics=_aggregate_diagnostics(point_diagnostics, model),
                    backend_identifier=spectrum.backend_identifier,
                    backend_version=spectrum.backend_version,
                )
            )
        return tuple(results)

    def evaluate_family(
        self,
        models: tuple[BlockHoppingModel, ...],
        fractional_wave_vectors: MatrixQuantity,
        tolerance: BlockHoppingBandTolerance,
    ) -> tuple[BlockHoppingBandEvaluation, ...]:
        """Evaluate compatible models with shared Bloch and spectrum batches.

        Parameters
        ----------
        models : tuple of BlockHoppingModel
            Nonempty ordered model family accepted by the Bloch evaluator.
        fractional_wave_vectors : MatrixQuantity
            Exact shared unitless reciprocal-coordinate record.
        tolerance : BlockHoppingBandTolerance
            Shared numerical Hermiticity tolerance.

        Returns
        -------
        tuple of BlockHoppingBandEvaluation
            One complete band result per model in supplied order.

        Raises
        ------
        TypeError
            If model, coordinate, tolerance, or delegated result types are invalid.
        ValueError
            If family metadata, result correlation, units, or Hermiticity gates fail.
        FloatingPointError
            If delegated Fourier or spectral arithmetic becomes nonfinite.
        RuntimeError
            If the complete Hermitian eigensolver does not converge.

        Notes
        -----
        The maintained default performs one shared-grid Fourier call and one
        concatenated complete-spectrum solve for a nonempty family and point grid.
        """
        if type(tolerance) is not BlockHoppingBandTolerance:
            raise TypeError("tolerance must be BlockHoppingBandTolerance")
        bloch_evaluations = self._bloch_evaluator.evaluate_family(
            models,
            fractional_wave_vectors,
        )
        if type(bloch_evaluations) is not tuple:
            raise TypeError("bloch_evaluator must return an exact tuple")
        if len(bloch_evaluations) != len(models):
            raise ValueError("Bloch family result count must match model count")
        for model, evaluation in zip(models, bloch_evaluations, strict=True):
            if type(evaluation) is not BlockHoppingBlochEvaluation:
                raise TypeError(
                    "bloch_evaluator family entries must be BlockHoppingBlochEvaluation"
                )
            if evaluation.model is not model:
                raise ValueError("Bloch family results must retain supplied models")
            if evaluation.fractional_wave_vectors is not fractional_wave_vectors:
                raise ValueError(
                    "Bloch family results must retain supplied coordinates"
                )
        return self.evaluate_bloch_family(bloch_evaluations, tolerance)

    def evaluate(
        self,
        model: BlockHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
        tolerance: BlockHoppingBandTolerance,
    ) -> BlockHoppingBandEvaluation:
        """The action evaluates Bloch matrices once and returns complete bands.

        Parameters
        ----------
        model : BlockHoppingModel
            Matrix-valued hopping model to evaluate.
        fractional_wave_vectors : MatrixQuantity
            Unitless fractional reciprocal coordinates with one row per point.
        tolerance : BlockHoppingBandTolerance
            Caller-supplied numerical Hermiticity tolerance.

        Returns
        -------
        BlockHoppingBandEvaluation
            Result retaining the exact Bloch evaluation produced by this action.

        Raises
        ------
        TypeError
            If the model, tolerance, coordinates, or delegated Bloch result has the
            wrong nominal type.
        ValueError
            If coordinate semantics, units, delegated input identities, or a
            nonempty-batch Hermiticity cutoff requirement is violated.
        FloatingPointError
            If delegated Bloch or band numerical arithmetic becomes nonfinite, or
            a constructed representative is not exactly Hermitian.
        RuntimeError
            If the complete Hermitian eigensolver does not converge.

        Notes
        -----
        This convenience action invokes :attr:`bloch_evaluator` exactly once and
        delegates the resulting object to :meth:`evaluate_bloch`.
        """
        if type(model) is not BlockHoppingModel:
            raise TypeError("model must be BlockHoppingModel")
        if type(tolerance) is not BlockHoppingBandTolerance:
            raise TypeError("tolerance must be BlockHoppingBandTolerance")
        if tolerance.hermiticity_absolute.unit != model.energy_unit:
            raise ValueError("absolute tolerance must use the model energy unit")
        bloch_evaluation = self._bloch_evaluator.evaluate(
            model,
            fractional_wave_vectors,
        )
        if type(bloch_evaluation) is not BlockHoppingBlochEvaluation:
            raise TypeError("bloch_evaluator must return BlockHoppingBlochEvaluation")
        if bloch_evaluation.model is not model:
            raise ValueError("Bloch evaluation must retain the supplied model")
        if bloch_evaluation.fractional_wave_vectors is not fractional_wave_vectors:
            raise ValueError("Bloch evaluation must retain the supplied coordinates")
        return self.evaluate_bloch(bloch_evaluation, tolerance)


class BlockHoppingRangeBandEvaluator:
    """Evaluate one source and an explicit ordered truncation family efficiently.

    Parameters
    ----------
    band_evaluator : BlockHoppingBandEvaluator or None, optional
        Compatible band evaluator. ``None`` constructs the maintained default.
    truncator : ConjugateClosedBlockHoppingTruncator or None, optional
        Policy-free ordinary-lattice truncator. ``None`` constructs the maintained
        default.

    Notes
    -----
    The action chooses no cutoff or range metric. It shares Fourier phases across the
    source and every retained model, then concatenates all nonempty Bloch matrices into
    one canonical complete-spectrum solve.
    """

    __slots__ = ("_band_evaluator", "_truncator")

    def __init__(
        self,
        band_evaluator: BlockHoppingBandEvaluator | None = None,
        truncator: ConjugateClosedBlockHoppingTruncator | None = None,
    ) -> None:
        """Retain injected nominal actions or construct maintained defaults."""
        if band_evaluator is None:
            band_evaluator = BlockHoppingBandEvaluator()
        elif not isinstance(band_evaluator, BlockHoppingBandEvaluator):
            raise TypeError("band_evaluator must be BlockHoppingBandEvaluator or None")
        if truncator is None:
            truncator = ConjugateClosedBlockHoppingTruncator()
        elif not isinstance(truncator, ConjugateClosedBlockHoppingTruncator):
            raise TypeError(
                "truncator must be ConjugateClosedBlockHoppingTruncator or None"
            )
        self._band_evaluator = band_evaluator
        self._truncator = truncator

    @property
    def band_evaluator(self) -> BlockHoppingBandEvaluator:
        """Return the retained complete-band evaluator."""
        return self._band_evaluator

    @property
    def truncator(self) -> ConjugateClosedBlockHoppingTruncator:
        """Return the retained conjugate-closed truncator."""
        return self._truncator

    def evaluate(
        self,
        source: BlockHoppingModel,
        selections: tuple[BlockHoppingRangeBandSelection, ...],
        fractional_wave_vectors: MatrixQuantity,
        tolerance: BlockHoppingBandTolerance,
    ) -> BlockHoppingRangeBandEvaluation:
        """Evaluate source and explicitly selected retained models as one family.

        Parameters
        ----------
        source : BlockHoppingModel
            Ordinary-integer-lattice source model retained by every truncation result.
        selections : tuple of BlockHoppingRangeBandSelection
            Nonempty ordered inventories and distinct retained-model identifiers.
        fractional_wave_vectors : MatrixQuantity
            Shared unitless reciprocal coordinates.
        tolerance : BlockHoppingBandTolerance
            Shared numerical Hermiticity tolerance.

        Returns
        -------
        BlockHoppingRangeBandEvaluation
            Source result and ordered exact truncation/band pairs.

        Raises
        ------
        TypeError
            If source, selections, coordinates, or tolerance have invalid nominal
            ownership.
        ValueError
            If inventories or identifiers repeat, a selection is incompatible with
            the source, coordinate semantics fail, or a Hermiticity gate rejects.
        FloatingPointError
            If delegated Fourier, diagnostic, representative, or spectral arithmetic
            becomes nonfinite.
        RuntimeError
            If the complete Hermitian eigensolver does not converge.

        Notes
        -----
        Selection order determines entry order only. It does not imply increasing
        physical range, improving accuracy, or scientific convergence.
        """
        if type(source) is not BlockHoppingModel:
            raise TypeError("source must be BlockHoppingModel")
        if type(selections) is not tuple or any(
            type(selection) is not BlockHoppingRangeBandSelection
            for selection in selections
        ):
            raise TypeError(
                "selections must be a tuple of BlockHoppingRangeBandSelection"
            )
        if not selections:
            raise ValueError("selections must be nonempty")
        identifiers = tuple(selection.truncated_identifier for selection in selections)
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("truncated model identifiers must be unique")
        inventories = tuple(
            tuple(
                displacement.components
                for displacement in selection.retained_displacements
            )
            for selection in selections
        )
        if len(set(inventories)) != len(inventories):
            raise ValueError("retained displacement inventories must be unique")

        # Truncation results are constructed once and retained verbatim in the final
        # entries, including their raw omitted-power evidence and source-term identity.
        truncations = tuple(
            self._truncator.execute(
                source,
                selection.retained_displacements,
                selection.truncated_identifier,
            )
            for selection in selections
        )
        if any(
            type(truncation) is not BlockHoppingTruncationResult
            for truncation in truncations
        ):
            raise TypeError("truncator must return BlockHoppingTruncationResult")
        for selection, truncation in zip(selections, truncations, strict=True):
            if truncation.source is not source:
                raise ValueError("truncation result must retain the supplied source")
            if truncation.retained_displacements != selection.retained_displacements:
                raise ValueError(
                    "truncation result must represent the selected inventory"
                )
            if truncation.truncated.identifier != selection.truncated_identifier:
                raise ValueError(
                    "truncation result must retain the selected identifier"
                )
        models = (source,) + tuple(truncation.truncated for truncation in truncations)
        evaluations = self._band_evaluator.evaluate_family(
            models,
            fractional_wave_vectors,
            tolerance,
        )
        if type(evaluations) is not tuple or len(evaluations) != len(models):
            raise ValueError("band evaluator must return one result per family model")
        if any(
            type(evaluation) is not BlockHoppingBandEvaluation
            for evaluation in evaluations
        ):
            raise TypeError(
                "band evaluator family entries must be BlockHoppingBandEvaluation"
            )
        return BlockHoppingRangeBandEvaluation(
            source_band_evaluation=evaluations[0],
            entries=tuple(
                BlockHoppingRangeBandEntry(truncation, evaluation)
                for truncation, evaluation in zip(
                    truncations,
                    evaluations[1:],
                    strict=True,
                )
            ),
        )
