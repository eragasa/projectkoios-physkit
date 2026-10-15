from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.lattice3d import DirectLattice3D
from projectkoios.physkit.periodic.unit_cell.base import AtomicBasis, UnitCell
from projectkoios.physkit.units.quantities import (
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
)


def test_composes_direct_lattice_scale_and_atomic_basis(
    fcc_direct_lattice: DirectLattice3D,
    silicon_basis: AtomicBasis,
) -> None:
    lattice_parameter = ScalarQuantity(5.43, PhysicalUnit("angstrom"))

    cell = UnitCell(
        direct_lattice=fcc_direct_lattice,
        lattice_parameter=lattice_parameter,
        atomic_basis=silicon_basis,
    )

    assert cell.direct_lattice is not fcc_direct_lattice
    np.testing.assert_array_equal(cell.A.magnitude, fcc_direct_lattice.A)
    assert cell.lattice_parameter is lattice_parameter
    assert cell.atomic_basis is silicon_basis


def test_exposes_dimensionless_and_physical_vectors_as_columns(
    silicon_basis: AtomicBasis,
) -> None:
    cell = UnitCell(
        direct_lattice=DirectLattice3D(
            a1=np.array([1.0, 0.0, 0.0]),
            a2=np.array([0.2, 2.0, 0.0]),
            a3=np.array([0.3, 0.4, 3.0]),
        ),
        lattice_parameter=ScalarQuantity(2.0, PhysicalUnit("angstrom")),
        atomic_basis=silicon_basis,
    )
    expected_a = np.array(
        (
            (1.0, 0.2, 0.3),
            (0.0, 2.0, 0.4),
            (0.0, 0.0, 3.0),
        )
    )

    np.testing.assert_array_equal(cell.A.magnitude, expected_a)
    np.testing.assert_array_equal(cell.H.magnitude, 2.0 * expected_a)
    np.testing.assert_array_equal(cell.a2.magnitude, expected_a[:, 1])
    np.testing.assert_array_equal(cell.h2.magnitude, 2.0 * expected_a[:, 1])
    assert isinstance(cell.A.unit, Unitless)
    assert cell.H.unit.expression == "angstrom"


def test__UnitCell__init__recovers_metric_from_normalized_direct_basis(
    silicon_basis: AtomicBasis,
) -> None:
    physical_a = 5.2
    physical_b = 6.1
    physical_c = 7.3
    requested_angles = np.array([75.0, 80.0, 65.0])
    lattice_parameter = ScalarQuantity(
        physical_a,
        PhysicalUnit("angstrom"),
    )
    normalized_lattice = DirectLattice3D.from_lattice_parameters(
        a=1.0,
        b=physical_b / physical_a,
        c=physical_c / physical_a,
        alpha_degrees=75.0,
        beta_degrees=80.0,
        gamma_degrees=65.0,
    )

    cell = UnitCell(
        direct_lattice=normalized_lattice,
        lattice_parameter=lattice_parameter,
        atomic_basis=silicon_basis,
    )

    normalized_lengths = np.linalg.norm(cell.A.magnitude, axis=0)
    physical_vectors = cell.H.magnitude
    physical_lengths = np.linalg.norm(physical_vectors, axis=0)
    angle_cosines = np.array(
        [
            np.dot(physical_vectors[:, 1], physical_vectors[:, 2])
            / (physical_lengths[1] * physical_lengths[2]),
            np.dot(physical_vectors[:, 0], physical_vectors[:, 2])
            / (physical_lengths[0] * physical_lengths[2]),
            np.dot(physical_vectors[:, 0], physical_vectors[:, 1])
            / (physical_lengths[0] * physical_lengths[1]),
        ]
    )
    physical_angles = np.degrees(np.arccos(np.clip(angle_cosines, -1.0, 1.0)))

    assert cell.lattice_parameter is lattice_parameter
    np.testing.assert_allclose(
        normalized_lengths,
        [1.0, physical_b / physical_a, physical_c / physical_a],
    )
    np.testing.assert_allclose(
        physical_lengths,
        [physical_a, physical_b, physical_c],
    )
    np.testing.assert_allclose(physical_angles, requested_angles)


def test_copies_and_freezes_direct_lattice_representation(
    fcc_direct_lattice: DirectLattice3D,
    silicon_basis: AtomicBasis,
) -> None:
    cell = UnitCell(
        direct_lattice=fcc_direct_lattice,
        lattice_parameter=ScalarQuantity(5.43, PhysicalUnit("angstrom")),
        atomic_basis=silicon_basis,
    )

    assert cell.direct_lattice is not fcc_direct_lattice
    assert cell.A.magnitude[0, 0] == 0.5
    assert fcc_direct_lattice.A.flags.writeable is False
    with pytest.raises(ValueError):
        fcc_direct_lattice.A[0, 0] = 99.0
    with pytest.raises(ValueError):
        cell.direct_lattice.A[0, 0] = 99.0
    with pytest.raises(ValueError):
        cell.H.magnitude[0, 0] = 99.0


def test_accepts_any_physical_length_unit(
    fcc_direct_lattice: DirectLattice3D,
    silicon_basis: AtomicBasis,
) -> None:
    cell = UnitCell(
        direct_lattice=fcc_direct_lattice,
        lattice_parameter=ScalarQuantity(10.26121286, PhysicalUnit("bohr")),
        atomic_basis=silicon_basis,
    )

    assert cell.lattice_parameter.unit.expression == "bohr"


def test_rejects_non_length_lattice_parameter(
    fcc_direct_lattice: DirectLattice3D,
    silicon_basis: AtomicBasis,
) -> None:
    with pytest.raises(ValueError, match="length dimensionality"):
        UnitCell(
            direct_lattice=fcc_direct_lattice,
            lattice_parameter=ScalarQuantity(
                5.43,
                PhysicalUnit("electron_volt"),
            ),
            atomic_basis=silicon_basis,
        )
