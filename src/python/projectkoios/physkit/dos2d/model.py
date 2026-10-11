"""Two-dimensional equal-weight band density-of-states model."""

from dataclasses import dataclass

from projectkoios.physkit.dos.base import (
    BaseUniformBandDensityOfStates,
)
from projectkoios.physkit.units import ScalarQuantity


@dataclass(frozen=True, slots=True, kw_only=True)
class UniformBandDensityOfStates2D(BaseUniformBandDensityOfStates):
    """Normalize unique equal-weight samples over a complete two-dimensional BZ."""

    brillouin_zone_measure: ScalarQuantity
    state_degeneracy: float = 1.0

    def __post_init__(self) -> None:
        """Check two-dimensional normalization parameters."""
        self._check_parameters()

    @property
    def dimension(self) -> int:
        """Return two reciprocal-space dimensions."""
        return 2
