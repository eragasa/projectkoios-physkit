r"""Complete unit-aware bands from matrix-valued Bloch Hamiltonians.

The public adapter evaluates :class:`BlockHoppingModel` Bloch matrices, applies a
caller-supplied numerical Hermiticity tolerance, constructs the explicit Hermitian
representative :math:`(H + H^\dagger)/2`, and returns every ascending eigenvalue at
every supplied fractional reciprocal coordinate.

The generic unit-free spectrum operation is delegated to the pinned
``projectkoios.numerics.linalg.hermitian_spectrum`` owner. No scientific tolerance,
gap metric, cutoff policy, acceptance decision, cyclic displacement interpretation,
or finite-mesh Nyquist rule is owned here.

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
from projectkoios.physkit.periodic.lattice.hopping import BlockHoppingModel
from projectkoios.physkit.solidstate.tb.bloch import (
    BlockHoppingBlochEvaluation,
    BlockHoppingBlochEvaluator,
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
        """The action evaluates complete bands from an existing Bloch evaluation.

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
            If an argument is not owned by its declared nominal PhysKit type.
        ValueError
            If the absolute tolerance unit differs from the model energy unit or
            any nonempty-batch Hermiticity defect exceeds its applied cutoff.
        FloatingPointError
            If diagnostic, symmetrization, representative, or eigensolver output
            arithmetic becomes nonfinite, or the representative is not exactly
            Hermitian.
        RuntimeError
            If the complete Hermitian eigensolver does not converge.

        Notes
        -----
        Empty Bloch evaluations return empty energies and ``diagnostics=None``
        without invoking the numerical eigensolver.
        """
        if type(bloch_evaluation) is not BlockHoppingBlochEvaluation:
            raise TypeError("bloch_evaluation must be BlockHoppingBlochEvaluation")
        if type(tolerance) is not BlockHoppingBandTolerance:
            raise TypeError("tolerance must be BlockHoppingBandTolerance")
        model = bloch_evaluation.model
        if tolerance.hermiticity_absolute.unit != model.energy_unit:
            raise ValueError("absolute tolerance must use the model energy unit")
        # The adapter strips energy units and array-backend ownership at the solver
        # boundary while preserving point, basis-row, and basis-column ordering.
        matrices = tuple(
            tuple(tuple(complex(entry) for entry in row) for row in matrix.magnitude)
            for matrix in bloch_evaluation.bloch_matrices
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
            # Empty batches declare their retained basis width explicitly because no
            # matrix exists from which the canonical owner could infer it.
            definition = BatchedHermitianSpectrumDefinition.empty(
                model.basis_size,
                numerical_tolerance,
            )
        spectrum = self._spectrum_solver.solve(definition)
        # Unit-free solver diagnostics regain the model energy unit point by point.
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
            for index, diagnostic in enumerate(spectrum.matrix_diagnostics)
        )
        # Aggregate owner diagnostics use the same restoration as point diagnostics;
        # absence remains meaningful for an empty reciprocal-point inventory.
        aggregate = spectrum.diagnostics
        diagnostics = (
            None
            if aggregate is None
            else BlockHoppingBandDiagnostics(
                maximum_input_hermiticity_absolute_defect=ScalarQuantity(
                    aggregate.maximum_input_hermiticity_absolute_defect,
                    model.energy_unit,
                ),
                maximum_representative_adjustment_absolute_defect=ScalarQuantity(
                    aggregate.maximum_representative_adjustment_absolute_defect,
                    model.energy_unit,
                ),
                maximum_input_matrix_absolute_entry=ScalarQuantity(
                    aggregate.maximum_input_matrix_absolute_entry,
                    model.energy_unit,
                ),
                maximum_applied_hermiticity_cutoff=ScalarQuantity(
                    aggregate.maximum_applied_hermiticity_cutoff,
                    model.energy_unit,
                ),
            )
        )
        # Eigenvalues regain their physical unit and remain correlated with the exact
        # Bloch evaluation rather than copied model or coordinate records.
        return BlockHoppingBandEvaluation(
            bloch_evaluation=bloch_evaluation,
            band_energies=MatrixQuantity(
                np.asarray(spectrum.eigenvalues, dtype=np.float64).reshape(
                    len(spectrum.eigenvalues),
                    model.basis_size,
                ),
                model.energy_unit,
            ),
            tolerance=tolerance,
            point_diagnostics=point_diagnostics,
            diagnostics=diagnostics,
            backend_identifier=spectrum.backend_identifier,
            backend_version=spectrum.backend_version,
        )

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
