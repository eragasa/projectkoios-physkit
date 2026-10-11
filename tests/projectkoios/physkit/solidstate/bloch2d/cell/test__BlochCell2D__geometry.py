"""Two-dimensional Bloch direct-cell geometry verification."""

import numpy as np

from projectkoios.physkit.solidstate.bloch2d.cell import BlochCell2D
from projectkoios.physkit.units import MatrixQuantity, UnitSystem


def test__geometry__returns_dual_basis_area_and_metric() -> None:
    basis = np.array([[4.0, 1.0], [0.0, 3.0]])
    cell = BlochCell2D(
        direct_basis=MatrixQuantity(
            magnitude=basis,
            unit=UnitSystem.METAL.length_unit,
        ),
        unit_system=UnitSystem.METAL,
    )

    np.testing.assert_allclose(
        basis.T @ cell.reciprocal_basis.magnitude,
        2.0 * np.pi * np.eye(2),
        atol=2e-15,
    )
    np.testing.assert_allclose(cell.metric_tensor.magnitude, basis.T @ basis)
    assert cell.area.magnitude == 12.0
