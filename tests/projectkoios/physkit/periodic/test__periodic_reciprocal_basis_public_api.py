"""Public façade tests for reciprocal-mode basis identities."""

from projectkoios.physkit import periodic
from projectkoios.physkit.periodic.basis import (
    ReciprocalModeBasis1D,
    ReciprocalModeBasis3D,
)


def test__periodic_facade__exports_reciprocal_mode_basis_identities() -> None:
    """The supported façade retains the exact defining-module classes."""
    assert periodic.ReciprocalModeBasis1D is ReciprocalModeBasis1D
    assert periodic.ReciprocalModeBasis3D is ReciprocalModeBasis3D
    assert "ReciprocalModeBasis1D" in periodic.__all__
    assert "ReciprocalModeBasis3D" in periodic.__all__
