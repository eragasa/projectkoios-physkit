"""Public package façade checks for lattice geometry."""

from __future__ import annotations

import numpy as np

import projectkoios.physkit.periodic.lattice as lattice
from projectkoios.physkit.periodic.lattice.metric import DirectLatticeMetric


class TestLatticePublicApi:
    """Keep the supported lattice composition surface explicit."""

    def test__exports_direct_lattice_metric_with_defining_module_identity(self) -> None:
        assert lattice.DirectLatticeMetric is DirectLatticeMetric
        assert lattice.__all__ == [
            "BravaisLattice",
            "BravaisLatticeKind",
            "DirectLattice",
            "DirectLattice1D",
            "DirectLattice2D",
            "DirectLattice3D",
            "DirectLatticeMetric",
            "FirstBrillouinZone",
            "FirstBrillouinZone1D",
            "FirstBrillouinZone2D",
            "FirstBrillouinZone3D",
            "Lattice",
            "Lattice3D",
            "ReciprocalLattice",
            "ReciprocalLattice1D",
            "ReciprocalLattice2D",
            "ReciprocalLattice3D",
            "WignerSeitzCell",
            "WignerSeitzCell1D",
            "WignerSeitzCell2D",
            "WignerSeitzCell3D",
        ]

    def test__composes_exported_two_dimensional_lattice_and_metric(self) -> None:
        direct_lattice = lattice.DirectLattice2D(
            a1=np.asarray([2.0, 0.0]),
            a2=np.asarray([0.5, 3.0]),
        )

        metric = lattice.DirectLatticeMetric(direct_lattice)

        assert metric.dimension == 2
        assert metric.direct_lattice is not direct_lattice
        np.testing.assert_allclose(metric.primitive_basis, direct_lattice.A)
