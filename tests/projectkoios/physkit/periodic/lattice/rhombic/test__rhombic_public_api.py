"""Public façade tests for the rhombic crystallography package."""

import projectkoios.physkit.periodic.lattice.rhombic as rhombic
from projectkoios.physkit.periodic.lattice.rhombic.path import (
    UnitAreaRhombicPath2D,
    UnitAreaRhombicPathPoint2D,
)


def test__rhombic_facade__retains_defining_module_identities() -> None:
    """The cohesive package exports only its reviewed owner contracts."""
    assert rhombic.UnitAreaRhombicPath2D is UnitAreaRhombicPath2D
    assert rhombic.UnitAreaRhombicPathPoint2D is UnitAreaRhombicPathPoint2D
    assert rhombic.__all__ == [
        "UnitAreaRhombicPath2D",
        "UnitAreaRhombicPathPoint2D",
    ]
