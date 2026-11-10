"""Tests for reciprocal evaluation of linear ``sp3s*`` channel models."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.solidstate.tb.bloch import DenseBlockHoppingBlochEvaluator
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
    HomopolarSp3sStarIntegralChannelEnum,
)
from projectkoios.physkit.solidstate.tb.slater_koster.linear_bloch import (
    HomopolarSp3sStarLinearBlochEvaluation,
    HomopolarSp3sStarLinearBlochEvaluator,
)
from projectkoios.physkit.solidstate.tb.slater_koster.linear_parameterization import (
    HomopolarSp3sStarLinearParameterizationAssembler,
)
from projectkoios.physkit.solidstate.tb.slater_koster.parameters import (
    HomopolarSp3sStarParameterGroup,
)
from projectkoios.physkit.solidstate.tb.slater_koster.sp3s_star import (
    HomopolarSp3sStarTwoCenterIntegrals,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    ScalarQuantity,
    UnitSystem,
    VectorQuantity,
)


class TestHomopolarSp3sStarLinearBlochEvaluator:
    """Verify nested order, reciprocal periodicity, and direct-model parity."""

    def test__execute__reconstructs_direct_bloch_matrices_in_nested_order(self) -> None:
        definition = self._definition()
        parameterization = HomopolarSp3sStarLinearParameterizationAssembler().execute(
            definition=definition
        )
        points = MatrixQuantity(
            magnitude=np.asarray(
                (
                    (0.0, 0.0, 0.0),
                    (0.17, -0.23, 0.31),
                    (1.17, -1.23, 2.31),
                )
            ),
            unit=UnitSystem.NONDIMENSIONAL.length_unit,
        )

        evaluation = HomopolarSp3sStarLinearBlochEvaluator().execute(
            parameterization=parameterization,
            fractional_wave_vectors=points,
        )
        direct = DenseBlockHoppingBlochEvaluator().evaluate(
            HomopolarSp3sStarCellGaugeModelAssembler()
            .execute(definition=definition)
            .model,
            points,
        )

        assert evaluation.parameterization is parameterization
        assert evaluation.fractional_wave_vectors is points
        assert len(evaluation.bloch_matrices) == 1
        assert len(evaluation.bloch_matrices[0]) == 7
        assert all(len(samples) == 3 for samples in evaluation.bloch_matrices[0])
        for sample_index, expected in enumerate(direct.bloch_matrices):
            reconstructed = np.zeros((10, 10), dtype=np.complex128)
            for scalar, channel_matrices in zip(
                definition.parameter_groups[0].integrals.ordered_quantities,
                evaluation.bloch_matrices[0],
                strict=True,
            ):
                reconstructed += (
                    scalar.magnitude * channel_matrices[sample_index].magnitude
                )
            np.testing.assert_allclose(
                reconstructed,
                expected.magnitude,
                rtol=2.0e-15,
                atol=2.0e-15,
            )
        for channel_matrices in evaluation.bloch_matrices[0]:
            np.testing.assert_allclose(
                channel_matrices[1].magnitude,
                channel_matrices[2].magnitude,
                rtol=0.0,
                atol=2.0e-14,
            )

    def test__matrix__uses_exact_group_channel_and_sample_keys(self) -> None:
        parameterization = HomopolarSp3sStarLinearParameterizationAssembler().execute(
            definition=self._definition()
        )
        points = MatrixQuantity(
            magnitude=np.asarray(((0.0, 0.0, 0.0),)),
            unit=UnitSystem.NONDIMENSIONAL.length_unit,
        )
        evaluation = HomopolarSp3sStarLinearBlochEvaluator().execute(
            parameterization=parameterization,
            fractional_wave_vectors=points,
        )

        matrix = evaluation.matrix(
            parameter_group_identifier="shared_group",
            channel=HomopolarSp3sStarIntegralChannelEnum.pp_sigma,
            sample_index=0,
        )

        assert type(matrix) is ComplexMatrixQuantity
        assert matrix is evaluation.bloch_matrices[0][2][0]
        assert matrix.unit == UnitSystem.NONDIMENSIONAL.energy_unit
        assert matrix.magnitude.flags.writeable is False
        with pytest.raises(TypeError):
            evaluation.matrix(
                parameter_group_identifier="shared_group",
                channel=HomopolarSp3sStarIntegralChannelEnum.pp_sigma,
                sample_index=True,
            )
        with pytest.raises(ValueError, match="outside"):
            evaluation.matrix(
                parameter_group_identifier="shared_group",
                channel=HomopolarSp3sStarIntegralChannelEnum.pp_sigma,
                sample_index=1,
            )

    def test__public_construction_and_action_require_named_arguments(self) -> None:
        parameterization = HomopolarSp3sStarLinearParameterizationAssembler().execute(
            definition=self._definition()
        )
        points = MatrixQuantity(
            magnitude=np.asarray(((0.0, 0.0, 0.0),)),
            unit=UnitSystem.NONDIMENSIONAL.length_unit,
        )
        evaluation = HomopolarSp3sStarLinearBlochEvaluator().execute(
            parameterization=parameterization,
            fractional_wave_vectors=points,
        )

        with pytest.raises(TypeError):
            HomopolarSp3sStarLinearBlochEvaluation(
                parameterization,
                points,
                evaluation.bloch_matrices,
            )
        with pytest.raises(TypeError):
            HomopolarSp3sStarLinearBlochEvaluator(DenseBlockHoppingBlochEvaluator())
        with pytest.raises(TypeError):
            HomopolarSp3sStarLinearBlochEvaluator().execute(
                parameterization,
                points,
            )
        with pytest.raises(ValueError, match="nondimensional"):
            HomopolarSp3sStarLinearBlochEvaluation(
                parameterization=parameterization,
                fractional_wave_vectors=MatrixQuantity(
                    magnitude=np.asarray(((0.0, 0.0, 0.0),)),
                    unit=UnitSystem.METAL.length_unit,
                ),
                bloch_matrices=evaluation.bloch_matrices,
            )
        with pytest.raises(ValueError, match="sample_count, 3"):
            HomopolarSp3sStarLinearBlochEvaluation(
                parameterization=parameterization,
                fractional_wave_vectors=MatrixQuantity(
                    magnitude=np.asarray(((0.0, 0.0),)),
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                ),
                bloch_matrices=evaluation.bloch_matrices,
            )

    @classmethod
    def _definition(cls) -> HomopolarSp3sStarCellGaugeModelDefinition:
        group = cls._parameter_group()
        return HomopolarSp3sStarCellGaugeModelDefinition(
            identifier="linear-bloch-sp3s-star",
            basis=cls._basis(),
            parameter_groups=(group,),
            bonds=(
                cls._bond(
                    displacement=(0, 0, 0),
                    cartesian_displacement=np.asarray((1.0, 1.0, 1.0)),
                ),
                cls._bond(
                    displacement=(1, -1, 0),
                    cartesian_displacement=np.asarray((1.0, -1.0, 0.5)),
                ),
            ),
            energy_reference="explicit-zero-reference",
        )

    @staticmethod
    def _basis() -> HomopolarSp3sStarSpatialBasis:
        centers = tuple(
            SlaterKosterBasisCenter(
                identifier=identifier,
                fractional_coordinate=VectorQuantity(
                    magnitude=np.asarray(coordinate),
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                ),
            )
            for identifier, coordinate in (
                ("origin_center", (0.0, 0.0, 0.0)),
                ("quarter_diagonal_center", (0.25, 0.25, 0.25)),
            )
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
                display_identifier="linear-bloch-two-center-sp3s-star",
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
            parameter_group_identifier="shared_group",
        )

    @staticmethod
    def _parameter_group() -> HomopolarSp3sStarParameterGroup:
        magnitudes = (1.0, -2.0, 3.0, -4.0, 5.0, 0.0, 0.0)
        unit = UnitSystem.METAL.energy_unit
        quantities = tuple(
            ScalarQuantity(magnitude=value, unit=unit) for value in magnitudes
        )
        return HomopolarSp3sStarParameterGroup(
            identifier="shared_group",
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
