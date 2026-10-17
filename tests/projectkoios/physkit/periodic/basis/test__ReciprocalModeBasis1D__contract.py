"""Contract tests for the canonical one-dimensional reciprocal-mode basis."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.basis import ReciprocalModeBasis1D
from projectkoios.physkit.periodic.lattice import (
    DirectLattice1D,
    ReciprocalLattice1D,
)
from projectkoios.physkit.periodic.reciprocal import (
    ReciprocalLattice1D as LegacyReciprocalLattice1D,
)


def _reciprocal_lattice() -> ReciprocalLattice1D:
    return ReciprocalLattice1D.from_direct_lattice(DirectLattice1D(2.0))


def test__init__retains_canonical_lattice_and_constructs_ordered_vectors() -> None:
    """Mode vectors come from the supported reciprocal-lattice identity."""
    reciprocal = _reciprocal_lattice()
    source_modes = np.asarray([-2, 0, 3], dtype=np.int64)

    basis = ReciprocalModeBasis1D(reciprocal, source_modes)

    assert basis.reciprocal_lattice is reciprocal
    np.testing.assert_array_equal(basis.modes, source_modes)
    np.testing.assert_allclose(basis.vectors, reciprocal.vector(source_modes))
    assert basis.size == 3
    assert basis.modes.flags.writeable is False
    assert basis.vectors.flags.writeable is False
    assert source_modes.flags.writeable is True


def test__from_mode_limit__constructs_symmetric_inventory() -> None:
    """The convenience constructor retains modes from minus to plus cutoff."""
    reciprocal = _reciprocal_lattice()

    basis = ReciprocalModeBasis1D.from_mode_limit(reciprocal, 2)

    np.testing.assert_array_equal(basis.modes, [-2, -1, 0, 1, 2])
    np.testing.assert_allclose(
        basis.shifted_wavevectors(0.25),
        reciprocal.vector(basis.modes) + 0.25,
    )


def test__init__rejects_legacy_reciprocal_lattice_identity() -> None:
    """The unexported legacy lattice is not a second nominal owner."""
    direct = DirectLattice1D(2.0)
    legacy = LegacyReciprocalLattice1D(direct)

    with pytest.raises(TypeError, match="must be ReciprocalLattice1D"):
        ReciprocalModeBasis1D(legacy, [-1, 0, 1])  # type: ignore[arg-type]


def test__init__rejects_empty_or_duplicate_mode_inventory() -> None:
    """A finite basis requires at least one unique integer mode."""
    reciprocal = _reciprocal_lattice()

    with pytest.raises(ValueError, match="cannot be empty"):
        ReciprocalModeBasis1D(reciprocal, [])
    with pytest.raises(ValueError, match="must be unique"):
        ReciprocalModeBasis1D(reciprocal, [0, 0])
