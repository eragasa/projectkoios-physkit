r"""Assemble the linear ordinary-cell basis of homopolar ``sp3s*`` models.

For an admitted cell-gauge definition with parameter groups :math:`g` and the stable
seven-channel order :math:`\alpha`, this module constructs one nondimensional hopping
model :math:`B_{g,\alpha}` for every explicit pair. The represented energy-valued model
then satisfies

.. math::

   H(V)=\sum_{g,\alpha}V_{g,\alpha}B_{g,\alpha}.

The module owns the reusable linear representation only. It does not own onsite terms,
parameter fitting, design-matrix flattening, residual construction, sample policy,
material interpretation, evidence, or scientific acceptance.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

from projectkoios.numerics.numpy.typing import Complex128Matrix
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    DenseBlockHoppingModel,
    SolidStateUnit,
)
from projectkoios.physkit.units import UnitSystem

from .assembly import (
    HomopolarSp3sStarCellGaugeModelDefinition,
    _accumulate_hopping,
    _dense_terms,
)
from .channels import (
    HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER,
    HomopolarSp3sStarIntegralChannelEnum,
    HomopolarSp3sStarTwoCenterChannelBasisConstructor,
)


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class HomopolarSp3sStarLinearParameterization(ResultsObject):
    r"""Correlate one structured definition with its channel hopping models.

    Parameters
    ----------
    definition
        Exact admitted cell-gauge model definition. Its parameter-group order is the
        outer order of ``channel_models``.
    channel_models
        Nested tuple in parameter-group then stable channel order. Every element is a
        nondimensional `DenseBlockHoppingModel` on the definition's exact basis.

    Raises
    ------
    TypeError
        If the definition, nested tuples, or models lack exact nominal ownership.
    ValueError
        If tuple lengths, model metadata, basis size, units, or deterministic model
        identifiers disagree with the retained definition.

    Notes
    -----
    The parameterization retains all seven channels for every declared group, including
    channels whose current energy coefficient is exactly zero. Zero values remain
    explicit coefficient choices and are never inferred from a group identifier.
    """

    definition: HomopolarSp3sStarCellGaugeModelDefinition
    channel_models: tuple[tuple[DenseBlockHoppingModel, ...], ...]

    def __post_init__(self) -> None:
        """Validate exact group/channel order and model metadata correlation."""
        if type(self.definition) is not HomopolarSp3sStarCellGaugeModelDefinition:
            raise TypeError(
                "definition must be HomopolarSp3sStarCellGaugeModelDefinition"
            )
        if type(self.channel_models) is not tuple:
            raise TypeError("channel_models must be a tuple")
        if len(self.channel_models) != len(self.definition.parameter_groups):
            raise ValueError("channel_models must match the parameter-group count")

        nondimensional_unit = UnitSystem.NONDIMENSIONAL.energy_unit
        expected_channel_count = len(HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER)
        for group, models in zip(
            self.definition.parameter_groups,
            self.channel_models,
            strict=True,
        ):
            if type(models) is not tuple:
                raise TypeError(
                    "each parameter-group channel inventory must be a tuple"
                )
            if len(models) != expected_channel_count:
                raise ValueError(
                    "each parameter group must contain seven channel models"
                )
            if any(type(model) is not DenseBlockHoppingModel for model in models):
                raise TypeError(
                    "channel_models must contain DenseBlockHoppingModel values"
                )
            for channel, model in zip(
                HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER,
                models,
                strict=True,
            ):
                if model.identifier != _linear_model_identifier(
                    definition_identifier=self.definition.identifier,
                    parameter_group_identifier=group.identifier,
                    channel=channel,
                ):
                    raise ValueError("channel-model identifier does not match its pair")
                if model.dimension is not LatticeDimension.THREE:
                    raise ValueError("channel-model dimension must be three")
                if model.energy_unit != nondimensional_unit:
                    raise ValueError("channel models must be nondimensional")
                if model.energy_reference != self.definition.energy_reference:
                    raise ValueError(
                        "channel-model energy_reference must match the definition"
                    )
                if (
                    model.basis_identifier
                    != self.definition.basis.definition.display_identifier
                ):
                    raise ValueError(
                        "channel-model basis_identifier must match the definition"
                    )
                if model.basis_size != self.definition.basis_size:
                    raise ValueError(
                        "channel-model basis_size must match the definition"
                    )

    @property
    def parameter_group_identifiers(self) -> tuple[str, ...]:
        """Return exact outer-order parameter-group identifiers."""
        return tuple(group.identifier for group in self.definition.parameter_groups)

    @property
    def channel_order(self) -> tuple[HomopolarSp3sStarIntegralChannelEnum, ...]:
        """Return the exact inner-order seven-channel enumeration."""
        return HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER

    def model(
        self,
        *,
        parameter_group_identifier: str,
        channel: HomopolarSp3sStarIntegralChannelEnum,
    ) -> DenseBlockHoppingModel:
        """Return one exactly identified nondimensional channel model.

        Parameters
        ----------
        parameter_group_identifier
            Built-in string matched by exact equality.
        channel
            Exact stable channel enum member.

        Returns
        -------
        DenseBlockHoppingModel
            Nondimensional ordinary-cell hopping model for the requested pair.

        Raises
        ------
        TypeError
            If either lookup key lacks its exact nominal owner.
        ValueError
            If the group identifier does not name exactly one declared group.
        """
        if type(parameter_group_identifier) is not str:
            raise TypeError("parameter_group_identifier must be a built-in string")
        if type(channel) is not HomopolarSp3sStarIntegralChannelEnum:
            raise TypeError("channel must be HomopolarSp3sStarIntegralChannelEnum")
        try:
            group_index = self.parameter_group_identifiers.index(
                parameter_group_identifier
            )
        except ValueError as error:
            raise ValueError(
                "parameter_group_identifier must name exactly one declared group"
            ) from error
        channel_index = HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER.index(channel)
        return self.channel_models[group_index][channel_index]


class HomopolarSp3sStarLinearParameterizationAssembler:
    r"""Assemble one nondimensional cell-gauge model per group and channel.

    Notes
    -----
    For each explicit pair :math:`(g,\alpha)`, the assembler selects only bonds naming
    :math:`g`, constructs :math:`B_\alpha(d_b)`, embeds the source-target block, and
    generates its Hermitian reverse exactly once. Parameter values are not read during
    this construction.
    """

    __slots__ = ()

    def execute(
        self,
        *,
        definition: HomopolarSp3sStarCellGaugeModelDefinition,
    ) -> HomopolarSp3sStarLinearParameterization:
        """Construct the exact nested linear channel-model inventory.

        Parameters
        ----------
        definition
            Validated structured model inventory. The argument is keyword-only.

        Returns
        -------
        HomopolarSp3sStarLinearParameterization
            Correlated definition and group-major, channel-minor model inventory.

        Raises
        ------
        TypeError
            If ``definition`` lacks its exact nominal owner.
        OverflowError
            If finite channel-block accumulation produces nonfinite values.
        """
        if type(definition) is not HomopolarSp3sStarCellGaugeModelDefinition:
            raise TypeError(
                "definition must be HomopolarSp3sStarCellGaugeModelDefinition"
            )

        channel_constructor = HomopolarSp3sStarTwoCenterChannelBasisConstructor()
        nondimensional_unit = cast(
            SolidStateUnit,
            UnitSystem.NONDIMENSIONAL.energy_unit,
        )
        grouped_models: list[tuple[DenseBlockHoppingModel, ...]] = []

        for group in definition.parameter_groups:
            models: list[DenseBlockHoppingModel] = []
            group_bonds = tuple(
                bond
                for bond in definition.bonds
                if bond.parameter_group_identifier == group.identifier
            )
            for channel in HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER:
                coefficients: dict[LatticeDisplacement, Complex128Matrix] = {}
                for parameterized_bond in group_bonds:
                    geometry = parameterized_bond.geometry
                    channel_matrix = channel_constructor.execute(
                        bond_direction=geometry.bond_direction
                    ).matrix(channel=channel)
                    _accumulate_hopping(
                        coefficients=coefficients,
                        definition=definition,
                        directed_bond=geometry.directed_bond,
                        hopping=channel_matrix.magnitude,
                    )
                models.append(
                    DenseBlockHoppingModel(
                        identifier=_linear_model_identifier(
                            definition_identifier=definition.identifier,
                            parameter_group_identifier=group.identifier,
                            channel=channel,
                        ),
                        dimension=LatticeDimension.THREE,
                        terms=_dense_terms(
                            coefficients=coefficients,
                            unit=nondimensional_unit,
                        ),
                        energy_unit=nondimensional_unit,
                        energy_reference=definition.energy_reference,
                        basis_identifier=definition.basis.definition.display_identifier,
                    )
                )
            grouped_models.append(tuple(models))

        return HomopolarSp3sStarLinearParameterization(
            definition=definition,
            channel_models=tuple(grouped_models),
        )


def _linear_model_identifier(
    *,
    definition_identifier: str,
    parameter_group_identifier: str,
    channel: HomopolarSp3sStarIntegralChannelEnum,
) -> str:
    """Return the deterministic identifier for one explicit linear pair."""
    return (
        f"{definition_identifier}:linear-channel:"
        f"{parameter_group_identifier}:{channel.value}"
    )
