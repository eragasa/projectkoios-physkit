"""Tests for homopolar ``sp3s*`` linear cell-model assembly."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.solidstate.tb.slater_koster.assembly import (
    HomopolarSp3sStarCellGaugeModelAssembler,
    HomopolarSp3sStarCellGaugeModelDefinition,
)
from projectkoios.physkit.solidstate.tb.slater_koster.basis import (
    HomopolarSp3sStarSpatialBasis,
    RealCubicOrbitalComponentEnum,
    SlaterKosterBasisCenter,
    SlaterKosterSpatialBasisCoordinate,
    SlaterKosterSpatialBasisDefinition,
)
from projectkoios.physkit.solidstate.tb.slater_koster.bond import (
    HomopolarSp3sStarBondGeometry,
    HomopolarSp3sStarParameterizedBond,
    SlaterKosterDirectedBond,
)
from projectkoios.physkit.solidstate.tb.slater_koster.channels import (
    HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER,
    HomopolarSp3sStarIntegralChannelEnum,
)
from projectkoios.physkit.solidstate.tb.slater_koster.linear_parameterization import (
    HomopolarSp3sStarLinearParameterization,
    HomopolarSp3sStarLinearParameterizationAssembler,
)
from projectkoios.physkit.solidstate.tb.slater_koster.parameters import (
    HomopolarSp3sStarParameterGroup,
)
from projectkoios.physkit.solidstate.tb.slater_koster.sp3s_star import (
    HomopolarSp3sStarTwoCenterIntegrals,
)
from projectkoios.physkit.units import ScalarQuantity, UnitSystem, VectorQuantity


class TestHomopolarSp3sStarLinearParameterizationAssembler:
    """Verify exact group/channel ordering and direct-model reconstruction."""

    def test__execute__reconstructs_direct_cell_model_with_explicit_zero_channels(
        self,
    ) -> None:
        definition = self._definition()

        parameterization = HomopolarSp3sStarLinearParameterizationAssembler().execute(
            definition=definition
        )
        direct = (
            HomopolarSp3sStarCellGaugeModelAssembler()
            .execute(definition=definition)
            .model
        )

        assert parameterization.definition is definition
        assert parameterization.parameter_group_identifiers == (
            "first_group",
            "second_group",
        )
        assert parameterization.channel_order == HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER
        reconstructed = self._reconstruct_coefficients(parameterization)
        assert tuple(sorted(reconstructed)) == tuple(
            term.displacement.components for term in direct.terms
        )
        for term in direct.terms:
            np.testing.assert_allclose(
                reconstructed[term.displacement.components],
                term.coefficient.magnitude,
                rtol=0.0,
                atol=2.0e-15,
            )
        # The first group deliberately fixes two channels to exact zero. Their basis
        # models remain present, while their scalar coefficients contribute nothing.
        first_integrals = definition.parameter_groups[0].integrals.ordered_quantities
        assert first_integrals[5].magnitude == 0.0
        assert first_integrals[6].magnitude == 0.0
        assert len(parameterization.channel_models[0]) == 7

    def test__execute__keeps_parameter_groups_independent(self) -> None:
        parameterization = HomopolarSp3sStarLinearParameterizationAssembler().execute(
            definition=self._definition()
        )

        first = parameterization.model(
            parameter_group_identifier="first_group",
            channel=HomopolarSp3sStarIntegralChannelEnum.ss_sigma,
        )
        second = parameterization.model(
            parameter_group_identifier="second_group",
            channel=HomopolarSp3sStarIntegralChannelEnum.ss_sigma,
        )

        assert tuple(term.displacement.components for term in first.terms) == (
            (0, 0, 0),
        )
        assert tuple(term.displacement.components for term in second.terms) == (
            (-1, 0, 0),
            (1, 0, 0),
        )
        assert first.energy_unit == UnitSystem.NONDIMENSIONAL.energy_unit
        assert second.energy_unit == UnitSystem.NONDIMENSIONAL.energy_unit
        assert all(
            term.coefficient.magnitude.flags.writeable is False
            for model in (first, second)
            for term in model.terms
        )

    def test__public_actions_and_lookup_reject_positional_or_inferred_keys(
        self,
    ) -> None:
        definition = self._definition()
        parameterization = HomopolarSp3sStarLinearParameterizationAssembler().execute(
            definition=definition
        )

        with pytest.raises(TypeError):
            HomopolarSp3sStarLinearParameterization(
                definition,
                parameterization.channel_models,
            )
        with pytest.raises(TypeError):
            HomopolarSp3sStarLinearParameterizationAssembler().execute(definition)
        with pytest.raises(TypeError):
            parameterization.model(
                parameter_group_identifier="first_group",
                channel="ss_sigma",
            )
        with pytest.raises(ValueError, match="exactly one declared group"):
            parameterization.model(
                parameter_group_identifier="undeclared_group",
                channel=HomopolarSp3sStarIntegralChannelEnum.ss_sigma,
            )

    @classmethod
    def _reconstruct_coefficients(
        cls,
        parameterization: HomopolarSp3sStarLinearParameterization,
    ) -> dict[tuple[int, int, int], np.ndarray]:
        coefficients: dict[tuple[int, int, int], np.ndarray] = {}
        for group, models in zip(
            parameterization.definition.parameter_groups,
            parameterization.channel_models,
            strict=True,
        ):
            for scalar, model in zip(
                group.integrals.ordered_quantities,
                models,
                strict=True,
            ):
                for term in model.terms:
                    displacement = term.displacement.components
                    if displacement not in coefficients:
                        coefficients[displacement] = np.zeros(
                            term.coefficient.magnitude.shape,
                            dtype=np.complex128,
                        )
                    coefficients[displacement] += (
                        scalar.magnitude * term.coefficient.magnitude
                    )
        return coefficients

    @classmethod
    def _definition(cls) -> HomopolarSp3sStarCellGaugeModelDefinition:
        return HomopolarSp3sStarCellGaugeModelDefinition(
            identifier="linear-two-group-sp3s-star",
            basis=cls._basis(),
            parameter_groups=(
                cls._parameter_group(
                    identifier="first_group",
                    magnitudes=(1.0, -2.0, 3.0, -4.0, 5.0, 0.0, 0.0),
                ),
                cls._parameter_group(
                    identifier="second_group",
                    magnitudes=(-0.5, 0.75, 1.25, -1.5, 2.0, -2.25, 2.5),
                ),
            ),
            bonds=(
                cls._bond(
                    displacement=(0, 0, 0),
                    cartesian_displacement=np.asarray((1.0, 1.0, 1.0)),
                    parameter_group_identifier="first_group",
                ),
                cls._bond(
                    displacement=(1, 0, 0),
                    cartesian_displacement=np.asarray((1.0, -1.0, 0.5)),
                    parameter_group_identifier="second_group",
                ),
            ),
            energy_reference="explicit-zero-reference",
        )

    @staticmethod
    def _basis() -> HomopolarSp3sStarSpatialBasis:
        centers = (
            SlaterKosterBasisCenter(
                identifier="origin_center",
                fractional_coordinate=VectorQuantity(
                    magnitude=np.asarray((0.0, 0.0, 0.0)),
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                ),
            ),
            SlaterKosterBasisCenter(
                identifier="quarter_diagonal_center",
                fractional_coordinate=VectorQuantity(
                    magnitude=np.asarray((0.25, 0.25, 0.25)),
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                ),
            ),
        )
        block = (
            (RealCubicOrbitalComponentEnum.s, 0),
            (RealCubicOrbitalComponentEnum.p_x, 0),
            (RealCubicOrbitalComponentEnum.p_y, 0),
            (RealCubicOrbitalComponentEnum.p_z, 0),
            (RealCubicOrbitalComponentEnum.s, 1),
        )
        return HomopolarSp3sStarSpatialBasis(
            definition=SlaterKosterSpatialBasisDefinition(
                display_identifier="linear-two-center-sp3s-star",
                centers=centers,
                coordinates=tuple(
                    SlaterKosterSpatialBasisCoordinate(
                        center_identifier=center.identifier,
                        orbital_component=component,
                        radial_channel=radial_channel,
                    )
                    for center in centers
                    for component, radial_channel in block
                ),
            )
        )

    @staticmethod
    def _bond(
        *,
        displacement: tuple[int, int, int],
        cartesian_displacement: np.ndarray,
        parameter_group_identifier: str,
    ) -> HomopolarSp3sStarParameterizedBond:
        return HomopolarSp3sStarParameterizedBond(
            geometry=HomopolarSp3sStarBondGeometry(
                directed_bond=SlaterKosterDirectedBond(
                    ordinary_cell_displacement=LatticeDisplacement(
                        dimension=LatticeDimension.THREE,
                        components=displacement,
                    ),
                    source_center_identifier="origin_center",
                    target_center_identifier="quarter_diagonal_center",
                ),
                cartesian_displacement=VectorQuantity(
                    magnitude=cartesian_displacement,
                    unit=UnitSystem.METAL.length_unit,
                ),
                unit_system=UnitSystem.METAL,
            ),
            parameter_group_identifier=parameter_group_identifier,
        )

    @staticmethod
    def _parameter_group(
        *,
        identifier: str,
        magnitudes: tuple[float, float, float, float, float, float, float],
    ) -> HomopolarSp3sStarParameterGroup:
        unit = UnitSystem.METAL.energy_unit
        quantities = tuple(
            ScalarQuantity(magnitude=value, unit=unit) for value in magnitudes
        )
        return HomopolarSp3sStarParameterGroup(
            identifier=identifier,
            integrals=HomopolarSp3sStarTwoCenterIntegrals(
                ss_sigma=quantities[0],
                sp_sigma=quantities[1],
                pp_sigma=quantities[2],
                pp_pi=quantities[3],
                s_star_p_sigma=quantities[4],
                s_s_star_sigma=quantities[5],
                s_star_s_star_sigma=quantities[6],
                unit_system=UnitSystem.METAL,
            ),
        )
