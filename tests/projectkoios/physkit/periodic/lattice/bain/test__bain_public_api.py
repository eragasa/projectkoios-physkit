"""Public façade tests for the Bain crystallography package."""

import projectkoios.physkit.periodic.lattice.bain as bain
from projectkoios.physkit.periodic.lattice.bain.bridge import (
    FccBainCoordinateBridge3D,
)
from projectkoios.physkit.periodic.lattice.bain.path import (
    DiamondDecoratedBainPath3D,
    DiamondDecoratedBainPathPoint3D,
    UnitVolumeBainPath3D,
    UnitVolumeBainPathPoint3D,
)


def test__bain_facade__retains_defining_module_identities() -> None:
    """The cohesive package exports only its reviewed owner contracts."""
    assert bain.DiamondDecoratedBainPath3D is DiamondDecoratedBainPath3D
    assert bain.DiamondDecoratedBainPathPoint3D is DiamondDecoratedBainPathPoint3D
    assert bain.FccBainCoordinateBridge3D is FccBainCoordinateBridge3D
    assert bain.UnitVolumeBainPath3D is UnitVolumeBainPath3D
    assert bain.UnitVolumeBainPathPoint3D is UnitVolumeBainPathPoint3D
    assert bain.__all__ == [
        "DiamondDecoratedBainPath3D",
        "DiamondDecoratedBainPathPoint3D",
        "FccBainCoordinateBridge3D",
        "UnitVolumeBainPath3D",
        "UnitVolumeBainPathPoint3D",
    ]
