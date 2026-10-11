"""One-dimensional equal-weight band density-of-states model."""

from dataclasses import dataclass

from projectkoios.physkit.dos.base import (
    BaseUniformBandDensityOfStates,
)
from projectkoios.physkit.units import ScalarQuantity


@dataclass(frozen=True, slots=True, kw_only=True)
class UniformBandDensityOfStates1D(BaseUniformBandDensityOfStates):
    """Normalize unique equal-weight samples over a complete one-dimensional BZ."""

    brillouin_zone_measure: ScalarQuantity
    state_degeneracy: float = 1.0

    def __post_init__(self) -> None:
        """Check one-dimensional normalization parameters."""
        self._check_parameters()

    @property
    def dimension(self) -> int:
        """Return one reciprocal-space dimension."""
        return 1
