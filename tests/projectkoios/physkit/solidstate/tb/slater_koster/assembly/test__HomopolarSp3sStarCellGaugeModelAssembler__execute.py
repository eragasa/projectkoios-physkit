"""Tests for ordinary-cell homopolar ``sp3s*`` model assembly."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.solidstate.tb.slater_koster.assembly import (
    HomopolarSp3sStarCellGaugeModelAssembler,
    HomopolarSp3sStarCellGaugeModelAssembly,
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
from projectkoios.physkit.solidstate.tb.slater_koster.parameters import (
    HomopolarSp3sStarParameterGroup,
)
from projectkoios.physkit.solidstate.tb.slater_koster.sp3s_star import (
    HomopolarSp3sStarTwoCenterHoppingConstructor,
    HomopolarSp3sStarTwoCenterHoppingDefinition,
    HomopolarSp3sStarTwoCenterIntegrals,
)
from projectkoios.physkit.units import ScalarQuantity, UnitSystem, VectorQuantity


class TestHomopolarSp3sStarCellGaugeModelAssembler:
    """Verify inventory admission, block placement, reversal, and API contracts."""

    def test__execute__assembles_zero_cell_forward_and_reverse_once(self) -> None:
        bond = self._bond(
            displacement=(0, 0, 0),
            source_center_identifier="origin_center",
            target_center_identifier="quarter_diagonal_center",
            cartesian_displacement=np.asarray((1.0, 1.0, 1.0)),
        )
        definition = self._definition(bonds=(bond,))

        assembly = HomopolarSp3sStarCellGaugeModelAssembler().execute(
            definition=definition
        )

        assert assembly.definition is definition
        assert len(assembly.model.terms) == 1
        term = assembly.model.terms[0]
        assert term.displacement.components == (0, 0, 0)
        assert assembly.model.basis_size == 10
        assert assembly.model.basis_identifier == "diamond-like-two-center-sp3s-star"
        coefficient = term.coefficient.magnitude
        expected_hopping = self._hopping(bond=bond)
        np.testing.assert_array_equal(coefficient[0:5, 5:10], expected_hopping)
        np.testing.assert_array_equal(coefficient[5:10, 0:5], expected_hopping.T)
        np.testing.assert_array_equal(coefficient, coefficient.conjugate().T)
        np.testing.assert_array_equal(coefficient[0:5, 0:5], np.zeros((5, 5)))
        np.testing.assert_array_equal(coefficient[5:10, 5:10], np.zeros((5, 5)))
        assert coefficient.flags.writeable is False

    def test__execute__sorts_nonzero_forward_and_reverse_displacements(self) -> None:
        bond = self._bond(
            displacement=(1, -2, 0),
            source_center_identifier="origin_center",
            target_center_identifier="quarter_diagonal_center",
            cartesian_displacement=np.asarray((2.0, -1.0, 2.0)),
        )

        assembly = HomopolarSp3sStarCellGaugeModelAssembler().execute(
            definition=self._definition(bonds=(bond,))
        )

        assert tuple(term.displacement.components for term in assembly.model.terms) == (
            (-1, 2, 0),
            (1, -2, 0),
        )
        by_displacement = {
            term.displacement.components: term.coefficient.magnitude
            for term in assembly.model.terms
        }
        expected_hopping = self._hopping(bond=bond)
        np.testing.assert_array_equal(
            by_displacement[(1, -2, 0)][0:5, 5:10],
            expected_hopping,
        )
        np.testing.assert_array_equal(
            by_displacement[(-1, 2, 0)][5:10, 0:5],
            expected_hopping.T,
        )

    def test__definition__rejects_duplicate_directed_bonds(self) -> None:
        first = self._bond(
            displacement=(1, 0, 0),
            source_center_identifier="origin_center",
            target_center_identifier="quarter_diagonal_center",
            cartesian_displacement=np.asarray((1.0, 0.0, 0.0)),
        )
        duplicate = self._bond(
            displacement=(1, 0, 0),
            source_center_identifier="origin_center",
            target_center_identifier="quarter_diagonal_center",
            cartesian_displacement=np.asarray((0.0, 1.0, 0.0)),
        )

        with pytest.raises(ValueError, match="directed bonds must be unique"):
            self._definition(bonds=(first, duplicate))

    def test__definition__rejects_explicit_hermitian_reverse(self) -> None:
        forward = self._bond(
            displacement=(1, 0, 0),
            source_center_identifier="origin_center",
            target_center_identifier="quarter_diagonal_center",
            cartesian_displacement=np.asarray((1.0, 0.0, 0.0)),
        )
        reverse = self._bond(
            displacement=(-1, 0, 0),
            source_center_identifier="quarter_diagonal_center",
            target_center_identifier="origin_center",
            cartesian_displacement=np.asarray((-1.0, 0.0, 0.0)),
        )

        with pytest.raises(ValueError, match="omit Hermitian reverses"):
            self._definition(bonds=(forward, reverse))

    def test__definition__rejects_endpoint_absent_from_basis(self) -> None:
        bond = self._bond(
            displacement=(0, 0, 0),
            source_center_identifier="origin_center",
            target_center_identifier="undeclared_center",
            cartesian_displacement=np.asarray((1.0, 0.0, 0.0)),
        )

        with pytest.raises(ValueError, match="endpoint must reference"):
            self._definition(bonds=(bond,))

    def test__definition__requires_every_parameter_group_to_be_used(self) -> None:
        bond = self._bond(
            displacement=(0, 0, 0),
            source_center_identifier="origin_center",
            target_center_identifier="quarter_diagonal_center",
            cartesian_displacement=np.asarray((1.0, 0.0, 0.0)),
        )
        groups = (
            self._parameter_group(identifier="used_group"),
            self._parameter_group(identifier="unused_group"),
        )
        bond = HomopolarSp3sStarParameterizedBond(
            geometry=bond.geometry,
            parameter_group_identifier="used_group",
        )

        with pytest.raises(ValueError, match="every parameter group must be used"):
            self._definition(bonds=(bond,), parameter_groups=groups)

    def test__definition__requires_one_exact_energy_unit(self) -> None:
        first = self._bond(
            displacement=(1, 0, 0),
            source_center_identifier="origin_center",
            target_center_identifier="quarter_diagonal_center",
            cartesian_displacement=np.asarray((1.0, 0.0, 0.0)),
            parameter_group_identifier="metal_group",
        )
        second = self._bond(
            displacement=(0, 1, 0),
            source_center_identifier="origin_center",
            target_center_identifier="quarter_diagonal_center",
            cartesian_displacement=np.asarray((0.0, 1.0, 0.0)),
            parameter_group_identifier="atomic_group",
        )
        groups = (
            self._parameter_group(identifier="metal_group"),
            self._parameter_group(
                identifier="atomic_group",
                unit_system=UnitSystem.ATOMIC,
            ),
        )

        with pytest.raises(ValueError, match="one exact energy unit"):
            self._definition(
                bonds=(first, second),
                parameter_groups=groups,
            )

    def test__public_actions_and_lookup_reject_positional_arguments(self) -> None:
        bond = self._bond(
            displacement=(0, 0, 0),
            source_center_identifier="origin_center",
            target_center_identifier="quarter_diagonal_center",
            cartesian_displacement=np.asarray((1.0, 0.0, 0.0)),
        )
        definition = self._definition(bonds=(bond,))
        assembly = HomopolarSp3sStarCellGaugeModelAssembler().execute(
            definition=definition
        )

        with pytest.raises(TypeError):
            HomopolarSp3sStarCellGaugeModelDefinition(
                definition.identifier,
                definition.basis,
                definition.parameter_groups,
                definition.bonds,
                definition.energy_reference,
            )
        with pytest.raises(TypeError):
            HomopolarSp3sStarCellGaugeModelAssembly(
                assembly.definition,
                assembly.model,
            )
        with pytest.raises(TypeError):
            HomopolarSp3sStarCellGaugeModelAssembler().execute(definition)
        with pytest.raises(TypeError):
            definition.parameter_group("default_group")

    @classmethod
    def _definition(
        cls,
        *,
        bonds: tuple[HomopolarSp3sStarParameterizedBond, ...],
        parameter_groups: tuple[HomopolarSp3sStarParameterGroup, ...] | None = None,
    ) -> HomopolarSp3sStarCellGaugeModelDefinition:
        return HomopolarSp3sStarCellGaugeModelDefinition(
            identifier="assembled-sp3s-star-cell-gauge",
            basis=cls._basis(),
            parameter_groups=parameter_groups
            or (cls._parameter_group(identifier="default_group"),),
            bonds=bonds,
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
                display_identifier="diamond-like-two-center-sp3s-star",
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

    @classmethod
    def _bond(
        cls,
        *,
        displacement: tuple[int, int, int],
        source_center_identifier: str,
        target_center_identifier: str,
        cartesian_displacement: np.ndarray,
        parameter_group_identifier: str = "default_group",
    ) -> HomopolarSp3sStarParameterizedBond:
        return HomopolarSp3sStarParameterizedBond(
            geometry=HomopolarSp3sStarBondGeometry(
                directed_bond=SlaterKosterDirectedBond(
                    ordinary_cell_displacement=LatticeDisplacement(
                        dimension=LatticeDimension.THREE,
                        components=displacement,
                    ),
                    source_center_identifier=source_center_identifier,
                    target_center_identifier=target_center_identifier,
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
        unit_system: UnitSystem = UnitSystem.METAL,
    ) -> HomopolarSp3sStarParameterGroup:
        unit = unit_system.energy_unit
        quantities = tuple(
            ScalarQuantity(magnitude=float(value), unit=unit) for value in range(1, 8)
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
                unit_system=unit_system,
            ),
        )

    @classmethod
    def _hopping(
        cls,
        *,
        bond: HomopolarSp3sStarParameterizedBond,
    ) -> np.ndarray:
        return (
            HomopolarSp3sStarTwoCenterHoppingConstructor()
            .execute(
                definition=HomopolarSp3sStarTwoCenterHoppingDefinition(
                    bond_direction=bond.geometry.bond_direction,
                    integrals=cls._parameter_group(
                        identifier=bond.parameter_group_identifier
                    ).integrals,
                )
            )
            .matrix.magnitude
        )
