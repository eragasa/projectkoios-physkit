from __future__ import annotations

import unittest
from collections.abc import Callable

import numpy as np

from projectkoios.physkit.periodic.lattice.bravais3d import (
    BravaisLattice,
    BravaisLatticeKind,
)
from projectkoios.physkit.periodic.lattice.lattice3d import DirectLattice3D

Rows = tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
]


class BravaisLatticeTest(unittest.TestCase):
    def test_retains_a_free_direct_lattice_by_identity(self) -> None:
        direct_lattice = DirectLattice3D(
            a1=np.asarray((1.0, 0.0, 0.0)),
            a2=np.asarray((0.0, 2.0, 0.0)),
            a3=np.asarray((0.0, 0.0, 3.0)),
        )

        lattice = BravaisLattice.free(direct_lattice)

        self.assertIs(lattice.kind, BravaisLatticeKind.free)
        self.assertIs(lattice.direct_lattice, direct_lattice)

    def test_constructs_documented_cubic_conventions(self) -> None:
        cases: tuple[
            tuple[Callable[[], BravaisLattice], BravaisLatticeKind, Rows], ...
        ] = (
            (
                BravaisLattice.cubic_primitive,
                BravaisLatticeKind.cubic_primitive,
                ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            ),
            (
                BravaisLattice.cubic_face_centered,
                BravaisLatticeKind.cubic_face_centered,
                ((-0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (-0.5, 0.5, 0.0)),
            ),
            (
                BravaisLattice.cubic_body_centered,
                BravaisLatticeKind.cubic_body_centered,
                ((0.5, 0.5, 0.5), (-0.5, 0.5, 0.5), (-0.5, -0.5, 0.5)),
            ),
            (
                BravaisLattice.cubic_body_centered_symmetric,
                BravaisLatticeKind.cubic_body_centered_symmetric,
                ((-0.5, 0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, -0.5)),
            ),
        )
        for factory, kind, expected in cases:
            with self.subTest(kind=kind):
                lattice = factory()
                self.assertIs(lattice.kind, kind)
                _assert_rows(self, lattice, expected)

    def test_constructs_documented_axial_conventions(self) -> None:
        cases: tuple[tuple[BravaisLattice, BravaisLatticeKind, Rows], ...] = (
            (
                BravaisLattice.hexagonal_trigonal_primitive(2.0),
                BravaisLatticeKind.hexagonal_trigonal_primitive,
                (
                    (1.0, 0.0, 0.0),
                    (-0.5, np.sqrt(3.0) / 2.0, 0.0),
                    (0.0, 0.0, 2.0),
                ),
            ),
            (
                BravaisLattice.trigonal_rhombohedral_axis_c(0.0),
                BravaisLatticeKind.trigonal_rhombohedral_axis_c,
                (
                    (np.sqrt(0.5), -np.sqrt(1.0 / 6.0), np.sqrt(1.0 / 3.0)),
                    (0.0, 2.0 * np.sqrt(1.0 / 6.0), np.sqrt(1.0 / 3.0)),
                    (-np.sqrt(0.5), -np.sqrt(1.0 / 6.0), np.sqrt(1.0 / 3.0)),
                ),
            ),
            (
                BravaisLattice.trigonal_rhombohedral_axis_111(0.0),
                BravaisLatticeKind.trigonal_rhombohedral_axis_111,
                (
                    (-1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0),
                    (2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0),
                    (2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0),
                ),
            ),
            (
                BravaisLattice.tetragonal_primitive(2.0),
                BravaisLatticeKind.tetragonal_primitive,
                ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 2.0)),
            ),
            (
                BravaisLattice.tetragonal_body_centered(2.0),
                BravaisLatticeKind.tetragonal_body_centered,
                ((0.5, -0.5, 1.0), (0.5, 0.5, 1.0), (-0.5, -0.5, 1.0)),
            ),
        )
        for lattice, kind, expected in cases:
            with self.subTest(kind=kind):
                self.assertIs(lattice.kind, kind)
                _assert_rows(self, lattice, expected)

    def test_constructs_documented_orthorhombic_conventions(self) -> None:
        cases: tuple[tuple[BravaisLattice, BravaisLatticeKind, Rows], ...] = (
            (
                BravaisLattice.orthorhombic_primitive(2.0, 3.0),
                BravaisLatticeKind.orthorhombic_primitive,
                ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)),
            ),
            (
                BravaisLattice.orthorhombic_base_centered(2.0, 3.0),
                BravaisLatticeKind.orthorhombic_base_centered,
                ((0.5, 1.0, 0.0), (-0.5, 1.0, 0.0), (0.0, 0.0, 3.0)),
            ),
            (
                BravaisLattice.orthorhombic_base_centered_alternate(2.0, 3.0),
                BravaisLatticeKind.orthorhombic_base_centered_alternate,
                ((0.5, -1.0, 0.0), (0.5, 1.0, 0.0), (0.0, 0.0, 3.0)),
            ),
            (
                BravaisLattice.orthorhombic_one_face_base_centered(2.0, 3.0),
                BravaisLatticeKind.orthorhombic_one_face_base_centered,
                ((1.0, 0.0, 0.0), (0.0, 1.0, -1.5), (0.0, 1.0, 1.5)),
            ),
            (
                BravaisLattice.orthorhombic_face_centered(2.0, 3.0),
                BravaisLatticeKind.orthorhombic_face_centered,
                ((0.5, 0.0, 1.5), (0.5, 1.0, 0.0), (0.0, 1.0, 1.5)),
            ),
            (
                BravaisLattice.orthorhombic_body_centered(2.0, 3.0),
                BravaisLatticeKind.orthorhombic_body_centered,
                ((0.5, 1.0, 1.5), (-0.5, 1.0, 1.5), (-0.5, -1.0, 1.5)),
            ),
        )
        for lattice, kind, expected in cases:
            with self.subTest(kind=kind):
                self.assertIs(lattice.kind, kind)
                _assert_rows(self, lattice, expected)

    def test_constructs_documented_monoclinic_and_triclinic_conventions(
        self,
    ) -> None:
        cases: tuple[tuple[BravaisLattice, BravaisLatticeKind, Rows], ...] = (
            (
                BravaisLattice.monoclinic_primitive_unique_c(2.0, 3.0, 0.0),
                BravaisLatticeKind.monoclinic_primitive_unique_c,
                ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)),
            ),
            (
                BravaisLattice.monoclinic_primitive_unique_b(2.0, 3.0, 0.0),
                BravaisLatticeKind.monoclinic_primitive_unique_b,
                ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)),
            ),
            (
                BravaisLattice.monoclinic_base_centered_unique_c(2.0, 3.0, 0.0),
                BravaisLatticeKind.monoclinic_base_centered_unique_c,
                ((0.5, 0.0, -1.5), (0.0, 2.0, 0.0), (0.5, 0.0, 1.5)),
            ),
            (
                BravaisLattice.monoclinic_base_centered_unique_b(2.0, 3.0, 0.0),
                BravaisLatticeKind.monoclinic_base_centered_unique_b,
                ((0.5, 1.0, 0.0), (-0.5, 1.0, 0.0), (0.0, 0.0, 3.0)),
            ),
            (
                BravaisLattice.triclinic(2.0, 3.0, 0.0, 0.0, 0.0),
                BravaisLatticeKind.triclinic,
                ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)),
            ),
        )
        for lattice, kind, expected in cases:
            with self.subTest(kind=kind):
                self.assertIs(lattice.kind, kind)
                _assert_rows(self, lattice, expected)

    def test_rejects_invalid_ratios_and_angle_cosines(self) -> None:
        with self.assertRaisesRegex(ValueError, "c_over_a"):
            BravaisLattice.tetragonal_primitive(0.0)
        with self.assertRaisesRegex(ValueError, "cos_gamma"):
            BravaisLattice.trigonal_rhombohedral_axis_c(-0.5)
        with self.assertRaisesRegex(ValueError, "positive volume"):
            BravaisLattice.triclinic(1.0, 1.0, 0.9, 0.9, -0.9)


def _assert_rows(
    case: unittest.TestCase,
    lattice: BravaisLattice,
    expected: Rows,
) -> None:
    np.testing.assert_allclose(lattice.direct_lattice.a1, expected[0], atol=1.0e-14)
    np.testing.assert_allclose(lattice.direct_lattice.a2, expected[1], atol=1.0e-14)
    np.testing.assert_allclose(lattice.direct_lattice.a3, expected[2], atol=1.0e-14)


if __name__ == "__main__":
    unittest.main()
