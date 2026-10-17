"""Canonical reciprocal-lattice tests for model-specific band visualizers."""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from projectkoios.physkit.periodic.lattice import (
    DirectLattice1D,
    ReciprocalLattice1D,
)
from projectkoios.physkit.solidstate.electronic.results import (
    ElectronicBandStructureResult,
)
from projectkoios.physkit.solidstate.phonons.results import PhononBandStructureResult
from projectkoios.physkit.solidstate.visualization.electronic_bands import (
    ElectronicBandStructureVisualizer,
)
from projectkoios.physkit.solidstate.visualization.phonon_bands import (
    PhononBandStructureVisualizer,
)


def test__electronic_visualizer__normalizes_by_canonical_basis_vector() -> None:
    """The older electronic visualizer accepts the supported lattice identity."""
    reciprocal = ReciprocalLattice1D(2.0)
    result = ElectronicBandStructureResult(
        k_points=np.asarray([-1.0, 0.0, 1.0]),
        energies=np.asarray([[1.0], [0.0], [1.0]]),
        eigenvectors=np.ones((3, 1, 1), dtype=np.complex128),
        reciprocal_vectors=np.asarray([0.0]),
    )
    figure, axes = plt.subplots()

    returned = ElectronicBandStructureVisualizer(result, reciprocal).plot(ax=axes)

    assert returned is axes
    np.testing.assert_array_equal(axes.lines[0].get_xdata(), [-0.5, 0.0, 0.5])
    plt.close(figure)


def test__phonon_visualizer__normalizes_by_canonical_basis_vector() -> None:
    """The older phonon visualizer accepts the supported lattice identity."""
    reciprocal = ReciprocalLattice1D(4.0)
    result = PhononBandStructureResult(
        q_points=np.asarray([-2.0, 0.0, 2.0]),
        angular_frequencies=np.asarray([[1.0], [0.0], [1.0]]),
        polarization_vectors=np.ones((3, 1, 1), dtype=np.complex128),
    )
    figure, axes = plt.subplots()

    returned = PhononBandStructureVisualizer(result, reciprocal).plot(ax=axes)

    assert returned is axes
    np.testing.assert_array_equal(axes.lines[0].get_xdata(), [-0.5, 0.0, 0.5])
    plt.close(figure)


def test__visualizers__reject_noncanonical_lattice_identity() -> None:
    """The model-specific visualizers enforce the reciprocal-lattice owner."""
    direct = DirectLattice1D(2.0)
    electronic = ElectronicBandStructureResult(
        k_points=np.asarray([0.0]),
        energies=np.asarray([[0.0]]),
        eigenvectors=np.ones((1, 1, 1), dtype=np.complex128),
        reciprocal_vectors=np.asarray([0.0]),
    )
    phonon = PhononBandStructureResult(
        q_points=np.asarray([0.0]),
        angular_frequencies=np.asarray([[0.0]]),
        polarization_vectors=np.ones((1, 1, 1), dtype=np.complex128),
    )

    with pytest.raises(TypeError, match="must be ReciprocalLattice1D"):
        ElectronicBandStructureVisualizer(electronic, direct)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="must be ReciprocalLattice1D"):
        PhononBandStructureVisualizer(phonon, direct)  # type: ignore[arg-type]
