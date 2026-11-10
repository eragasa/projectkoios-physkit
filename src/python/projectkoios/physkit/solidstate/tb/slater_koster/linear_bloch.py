r"""Evaluate homopolar ``sp3s*`` linear channel models on reciprocal samples.

For ordered reciprocal coordinates :math:`q_j`, this module evaluates every
nondimensional channel model and retains matrices in parameter-group, channel, then
sample order:

.. math::

   H_{\mathrm{cell}}(q_j;V)
   =\sum_{g,\alpha}V_{g,\alpha}B_{g,\alpha}(q_j).

The result deliberately does not flatten matrix entries, construct a fitting residual,
select samples, solve for parameters, or establish model adequacy.
"""

from __future__ import annotations

from dataclasses import dataclass

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.solidstate.tb.bloch import DenseBlockHoppingBlochEvaluator
from projectkoios.physkit.units import ComplexMatrixQuantity, MatrixQuantity, UnitSystem

from .channels import (
    HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER,
    HomopolarSp3sStarIntegralChannelEnum,
)
from .linear_parameterization import HomopolarSp3sStarLinearParameterization


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class HomopolarSp3sStarLinearBlochEvaluation(ResultsObject):
    r"""Retain channel matrices in group, channel, and reciprocal-sample order.

    Parameters
    ----------
    parameterization
        Exact linear ordinary-cell parameterization that owns the channel models.
    fractional_wave_vectors
        Unitless reciprocal coordinates retained by every model evaluation.
    bloch_matrices
        Nested tuple with order ``[parameter_group][channel][sample]``. Every matrix is
        nondimensional and acts on the parameterization's exact basis.

    Raises
    ------
    TypeError
        If values or nested tuple entries lack exact nominal ownership.
    ValueError
        If group, channel, sample, shape, or unit correlation is inconsistent.
    """

    parameterization: HomopolarSp3sStarLinearParameterization
    fractional_wave_vectors: MatrixQuantity
    bloch_matrices: tuple[
        tuple[tuple[ComplexMatrixQuantity, ...], ...],
        ...,
    ]

    def __post_init__(self) -> None:
        """Validate exact nested order, sample count, matrix shape, and unit."""
        if type(self.parameterization) is not HomopolarSp3sStarLinearParameterization:
            raise TypeError(
                "parameterization must be HomopolarSp3sStarLinearParameterization"
            )
        if type(self.fractional_wave_vectors) is not MatrixQuantity:
            raise TypeError("fractional_wave_vectors must be MatrixQuantity")
        if self.fractional_wave_vectors.unit != UnitSystem.NONDIMENSIONAL.length_unit:
            raise ValueError("fractional_wave_vectors must be nondimensional")
        if self.fractional_wave_vectors.magnitude.shape[1] != 3:
            raise ValueError(
                "fractional_wave_vectors must have shape (sample_count, 3)"
            )
        if type(self.bloch_matrices) is not tuple:
            raise TypeError("bloch_matrices must be a tuple")
        if len(self.bloch_matrices) != len(
            self.parameterization.definition.parameter_groups
        ):
            raise ValueError("bloch_matrices must match the parameter-group count")

        channel_count = len(HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER)
        sample_count = self.fractional_wave_vectors.magnitude.shape[0]
        basis_size = self.parameterization.definition.basis_size
        nondimensional_unit = UnitSystem.NONDIMENSIONAL.energy_unit
        for group_matrices in self.bloch_matrices:
            if type(group_matrices) is not tuple:
                raise TypeError("each parameter-group Bloch inventory must be a tuple")
            if len(group_matrices) != channel_count:
                raise ValueError("each parameter group must contain seven channels")
            for channel_matrices in group_matrices:
                if type(channel_matrices) is not tuple or any(
                    type(matrix) is not ComplexMatrixQuantity
                    for matrix in channel_matrices
                ):
                    raise TypeError(
                        "channel samples must be tuples of ComplexMatrixQuantity values"
                    )
                if len(channel_matrices) != sample_count:
                    raise ValueError(
                        "every channel must match the reciprocal sample count"
                    )
                if any(
                    matrix.magnitude.shape != (basis_size, basis_size)
                    for matrix in channel_matrices
                ):
                    raise ValueError(
                        "Bloch matrix shapes must match the parameterization basis"
                    )
                if any(
                    matrix.unit != nondimensional_unit for matrix in channel_matrices
                ):
                    raise ValueError("linear Bloch matrices must be nondimensional")

    def matrix(
        self,
        *,
        parameter_group_identifier: str,
        channel: HomopolarSp3sStarIntegralChannelEnum,
        sample_index: int,
    ) -> ComplexMatrixQuantity:
        """Return one exactly indexed nondimensional Bloch matrix.

        Parameters
        ----------
        parameter_group_identifier
            Built-in string matched by exact equality.
        channel
            Exact stable channel enum member.
        sample_index
            Built-in integer index into the retained reciprocal-coordinate rows.

        Returns
        -------
        ComplexMatrixQuantity
            Immutable matrix for the requested group, channel, and sample.

        Raises
        ------
        TypeError
            If a lookup key lacks its exact nominal owner.
        ValueError
            If the group is undeclared or the sample index is outside the retained
            sample range.
        """
        if type(sample_index) is not int:
            raise TypeError("sample_index must be a built-in int")
        self.parameterization.model(
            parameter_group_identifier=parameter_group_identifier,
            channel=channel,
        )
        group_index = self.parameterization.parameter_group_identifiers.index(
            parameter_group_identifier
        )
        channel_index = HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER.index(channel)
        if sample_index < 0 or sample_index >= len(
            self.bloch_matrices[group_index][channel_index]
        ):
            raise ValueError("sample_index is outside the retained sample range")
        return self.bloch_matrices[group_index][channel_index][sample_index]


class HomopolarSp3sStarLinearBlochEvaluator:
    r"""Evaluate all linear channel models through PhysKit's dense Bloch owner.

    Parameters
    ----------
    evaluator
        Optional nominal dense block-hopping Bloch evaluator. When omitted, PhysKit's
        default shared-grid evaluator is constructed.
    """

    __slots__ = ("_evaluator",)

    def __init__(
        self,
        *,
        evaluator: DenseBlockHoppingBlochEvaluator | None = None,
    ) -> None:
        """Retain an injected evaluator or construct the default owner."""
        if evaluator is None:
            evaluator = DenseBlockHoppingBlochEvaluator()
        elif type(evaluator) is not DenseBlockHoppingBlochEvaluator:
            raise TypeError("evaluator must be DenseBlockHoppingBlochEvaluator or None")
        self._evaluator = evaluator

    @property
    def evaluator(self) -> DenseBlockHoppingBlochEvaluator:
        """Return the exact dense Bloch evaluator used by this action."""
        return self._evaluator

    def execute(
        self,
        *,
        parameterization: HomopolarSp3sStarLinearParameterization,
        fractional_wave_vectors: MatrixQuantity,
    ) -> HomopolarSp3sStarLinearBlochEvaluation:
        """Evaluate the complete group-major, channel-minor model family.

        Parameters
        ----------
        parameterization
            Exact linear channel-model inventory.
        fractional_wave_vectors
            Ordered unitless reciprocal coordinates. Arguments are keyword-only.

        Returns
        -------
        HomopolarSp3sStarLinearBlochEvaluation
            Correlated nested matrix inventory in group, channel, and sample order.

        Raises
        ------
        TypeError
            If either argument lacks its exact nominal owner.
        ValueError
            If reciprocal coordinates are dimensional or have incompatible width.
        FloatingPointError
            If delegated finite-Fourier arithmetic produces nonfinite values.
        """
        if type(parameterization) is not HomopolarSp3sStarLinearParameterization:
            raise TypeError(
                "parameterization must be HomopolarSp3sStarLinearParameterization"
            )
        if type(fractional_wave_vectors) is not MatrixQuantity:
            raise TypeError("fractional_wave_vectors must be MatrixQuantity")

        flat_models = tuple(
            model
            for group_models in parameterization.channel_models
            for model in group_models
        )
        flat_evaluations = self._evaluator.evaluate_family(
            flat_models,
            fractional_wave_vectors,
        )
        channel_count = len(HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER)
        grouped_matrices = tuple(
            tuple(
                flat_evaluations[
                    group_index * channel_count + channel_index
                ].bloch_matrices
                for channel_index in range(channel_count)
            )
            for group_index in range(len(parameterization.channel_models))
        )
        return HomopolarSp3sStarLinearBlochEvaluation(
            parameterization=parameterization,
            fractional_wave_vectors=fractional_wave_vectors,
            bloch_matrices=grouped_matrices,
        )
