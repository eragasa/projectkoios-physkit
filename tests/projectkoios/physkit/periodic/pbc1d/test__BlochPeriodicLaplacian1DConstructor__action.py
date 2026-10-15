"""Numerical verification for one-dimensional Bloch-periodic seams."""

import numpy as np

from projectkoios.physkit.periodic.pbc1d.grid import PeriodicFiniteDifferenceGrid1D
from projectkoios.physkit.periodic.pbc1d.laplacian import (
    BlochPeriodicLaplacian1D,
    BlochPeriodicLaplacian1DConstructionRequest,
    BlochPeriodicLaplacian1DConstructor,
)
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity


class TestBlochPeriodicLaplacian1DConstructor:
    """Verify seam phases, Hermiticity, and represented eigenvalues."""

    @staticmethod
    def _construct(
        *, point_count: int, length: float, wave_number: float
    ) -> BlochPeriodicLaplacian1D:
        grid = PeriodicFiniteDifferenceGrid1D(
            point_count=point_count,
            oriented_cell_length=ScalarQuantity(
                length,
                PhysicalUnit("angstrom"),
            ),
        )
        return BlochPeriodicLaplacian1DConstructor().action(
            request=BlochPeriodicLaplacian1DConstructionRequest(
                grid=grid,
                bloch_wave_number=ScalarQuantity(
                    wave_number,
                    PhysicalUnit("1 / angstrom"),
                ),
            )
        )

    def test_zero_twist_closes_the_periodic_stencil(self) -> None:
        result = self._construct(point_count=8, length=4.0, wave_number=0.0)
        represented = result.represented_laplacian.to_csr().toarray()
        inverse_spacing_squared = 1.0 / result.request.grid.spacing.magnitude**2

        assert represented[0, -1] == inverse_spacing_squared
        assert represented[-1, 0] == inverse_spacing_squared
        np.testing.assert_allclose(represented, represented.conjugate().T)

    def test_nonzero_twist_places_conjugate_seam_phases(self) -> None:
        wave_number = 0.31
        result = self._construct(
            point_count=9,
            length=5.2,
            wave_number=wave_number,
        )
        represented = result.represented_laplacian.to_csr().toarray()
        scale = 1.0 / result.request.grid.spacing.magnitude**2
        phase = np.exp(1.0j * wave_number * 5.2)

        np.testing.assert_allclose(represented[-1, 0], phase * scale)
        np.testing.assert_allclose(represented[0, -1], np.conjugate(phase) * scale)
        np.testing.assert_allclose(represented, represented.conjugate().T)

    def test_eigenvalues_match_the_centered_difference_symbol(self) -> None:
        point_count = 12
        length = 5.4
        wave_number = 0.17
        result = self._construct(
            point_count=point_count,
            length=length,
            wave_number=wave_number,
        )
        represented = result.represented_laplacian.to_csr().toarray()
        observed = np.linalg.eigvalsh(represented)
        spacing = length / point_count
        integer_modes = np.fft.fftfreq(point_count, d=1.0 / point_count)
        physical_wave_numbers = wave_number + 2.0 * np.pi * integer_modes / length
        expected = (
            -4.0 * np.sin(0.5 * physical_wave_numbers * spacing) ** 2 / spacing**2
        )

        np.testing.assert_allclose(
            np.sort(observed),
            np.sort(expected),
            rtol=2e-14,
            atol=2e-14,
        )

    def test_grid_uses_endpoint_excluded_positions(self) -> None:
        result = self._construct(point_count=5, length=4.0, wave_number=0.0)

        np.testing.assert_allclose(
            result.request.grid.positions.magnitude,
            np.array([0.0, 0.8, 1.6, 2.4, 3.2]),
        )
