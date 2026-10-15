"""Nondegenerate abrupt p-n junction model in explicit SI units."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from projectkoios.numerics.numpy.typing import (
    Float64Array as FloatArray,
)
from projectkoios.numerics.numpy.typing import (
    NumpyArrayLike,
)
from projectkoios.physkit.constants import SI


def _positive_scalar(value: float, *, name: str) -> None:
    if type(value) not in (int, float):
        raise TypeError(f"{name} must be a built-in real number")
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")


@dataclass(frozen=True, slots=True, kw_only=True)
class SemiconductorMaterial:
    """Store bulk nondegenerate-semiconductor parameters in SI units.

    Band energies are measured in electron-volts, effective densities of
    states in inverse cubic meters, mobilities in square meters per volt
    second, and temperature in kelvin.
    """

    band_gap: float
    electron_affinity: float
    relative_permittivity: float
    conduction_density_of_states: float
    valence_density_of_states: float
    electron_mobility: float
    hole_mobility: float
    temperature: float

    def __post_init__(self) -> None:
        """Check every material parameter."""
        for name in (
            "band_gap",
            "electron_affinity",
            "relative_permittivity",
            "conduction_density_of_states",
            "valence_density_of_states",
            "electron_mobility",
            "hole_mobility",
            "temperature",
        ):
            _positive_scalar(getattr(self, name), name=name)

    @property
    def permittivity(self) -> float:
        """Return absolute permittivity in farads per meter."""
        return self.relative_permittivity * SI.eps0

    @property
    def thermal_energy(self) -> float:
        """Return ``k_B T`` in electron-volts."""
        return SI.k_B * self.temperature / SI.q

    @property
    def intrinsic_carrier_density(self) -> float:
        """Return the bulk nondegenerate intrinsic density in m⁻³."""
        density_prefactor = np.sqrt(
            self.conduction_density_of_states * self.valence_density_of_states
        )
        return float(
            density_prefactor * np.exp(-self.band_gap / (2.0 * self.thermal_energy))
        )

    @property
    def intrinsic_level_from_valence_band(self) -> float:
        """Return ``E_i - E_V`` in electron-volts."""
        return float(
            self.thermal_energy
            * np.log(self.valence_density_of_states / self.intrinsic_carrier_density)
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class AbruptDopingProfile:
    """Store uniform acceptor and donor densities in inverse cubic meters."""

    acceptor_density: float
    donor_density: float

    def __post_init__(self) -> None:
        """Check both doping densities."""
        _positive_scalar(self.acceptor_density, name="acceptor_density")
        _positive_scalar(self.donor_density, name="donor_density")


@dataclass(frozen=True, slots=True, kw_only=True)
class DepletionWidths:
    """Store total, p-side, and n-side depletion widths in meters."""

    total: float
    p_side: float
    n_side: float


@dataclass(frozen=True, slots=True, kw_only=True)
class EnergyBands:
    """Store equilibrium band energies in electron-volts."""

    conduction: FloatArray
    valence: FloatArray
    intrinsic: FloatArray
    fermi: FloatArray


@dataclass(frozen=True, slots=True, kw_only=True)
class CarrierProfiles:
    """Store electron and hole number densities in inverse cubic meters."""

    electron: FloatArray
    hole: FloatArray


@dataclass(frozen=True, slots=True, kw_only=True)
class CurrentProfiles:
    """Store drift, diffusion, and total current densities in A/m²."""

    electron_drift: FloatArray
    electron_diffusion: FloatArray
    electron_total: FloatArray
    hole_drift: FloatArray
    hole_diffusion: FloatArray
    hole_total: FloatArray


@dataclass(frozen=True, slots=True, kw_only=True)
class PNJunctionEquilibrium:
    """Store one finite-difference depletion-approximation evaluation."""

    position: FloatArray
    charge_density: FloatArray
    electrostatic_potential: FloatArray
    electric_field: FloatArray
    bands: EnergyBands
    carriers: CarrierProfiles
    currents: CurrentProfiles


@dataclass(frozen=True, slots=True, kw_only=True)
class AbruptPNJunction:
    """Evaluate a one-dimensional abrupt p-n junction in SI units.

    The implementation preserves the retired notebook's depletion charge,
    Dirichlet-potential, nondegenerate-carrier, and drift-diffusion equations.
    It is a pedagogical reduced model rather than a self-consistent device
    simulator.
    """

    material: SemiconductorMaterial
    doping: AbruptDopingProfile

    def __post_init__(self) -> None:
        """Check material and doping ownership."""
        if not isinstance(self.material, SemiconductorMaterial):
            raise TypeError("material must be SemiconductorMaterial")
        if not isinstance(self.doping, AbruptDopingProfile):
            raise TypeError("doping must be AbruptDopingProfile")

    def built_in_potential(self) -> float:
        """Return the equilibrium built-in potential in volts."""
        intrinsic_density = self.material.intrinsic_carrier_density
        argument = (
            self.doping.acceptor_density
            * self.doping.donor_density
            / intrinsic_density**2
        )
        if argument <= 1.0:
            raise ValueError("doping product must exceed the intrinsic-density square")
        return float(self.material.thermal_energy * np.log(argument))

    def depletion_widths(
        self,
        *,
        applied_potential: float = 0.0,
    ) -> DepletionWidths:
        """Return depletion widths for the selected applied potential."""
        if type(applied_potential) not in (int, float):
            raise TypeError("applied_potential must be a built-in real number")
        if not math.isfinite(applied_potential):
            raise ValueError("applied_potential must be finite")
        effective_potential = self.built_in_potential() - applied_potential
        if effective_potential <= 0.0:
            raise ValueError("applied_potential must be below the built-in potential")
        acceptors = self.doping.acceptor_density
        donors = self.doping.donor_density
        total = np.sqrt(
            2.0
            * self.material.permittivity
            / SI.q
            * (acceptors + donors)
            / (acceptors * donors)
            * effective_potential
        )
        n_side = total * acceptors / (acceptors + donors)
        p_side = total * donors / (acceptors + donors)
        return DepletionWidths(
            total=float(total),
            p_side=float(p_side),
            n_side=float(n_side),
        )

    def position_grid(
        self,
        *,
        applied_potential: float = 0.0,
        p_side_factor: float = 1.5,
        n_side_factor: float = 1.5,
        point_count: int = 600,
    ) -> FloatArray:
        """Return a grid extending beyond both depletion edges."""
        _positive_scalar(p_side_factor, name="p_side_factor")
        _positive_scalar(n_side_factor, name="n_side_factor")
        if type(point_count) is not int:
            raise TypeError("point_count must be a built-in integer")
        if point_count < 3:
            raise ValueError("point_count must be at least three")
        widths = self.depletion_widths(applied_potential=applied_potential)
        return np.linspace(
            -p_side_factor * widths.p_side,
            n_side_factor * widths.n_side,
            point_count,
        )

    def charge_density(
        self,
        position: NumpyArrayLike,
        *,
        applied_potential: float = 0.0,
    ) -> FloatArray:
        """Return depletion-approximation space charge in C/m³."""
        positions = self._position_array(position)
        widths = self.depletion_widths(applied_potential=applied_potential)
        density = np.zeros_like(positions)
        p_mask = (positions >= -widths.p_side) & (positions < 0.0)
        n_mask = (positions > 0.0) & (positions <= widths.n_side)
        density[p_mask] = -SI.q * self.doping.acceptor_density
        density[n_mask] = SI.q * self.doping.donor_density
        return density

    def solve_equilibrium(
        self,
        *,
        applied_potential: float = 0.0,
        p_side_factor: float = 1.5,
        n_side_factor: float = 1.5,
        point_count: int = 600,
    ) -> PNJunctionEquilibrium:
        """Evaluate potential, bands, carriers, and current-density terms."""
        position = self.position_grid(
            applied_potential=applied_potential,
            p_side_factor=p_side_factor,
            n_side_factor=n_side_factor,
            point_count=point_count,
        )
        charge_density = self.charge_density(
            position,
            applied_potential=applied_potential,
        )
        potential, field = self._solve_poisson(
            position,
            charge_density,
            applied_potential=applied_potential,
        )
        bands = self._energy_bands(position, potential)
        carriers = self._carrier_profiles(bands)
        currents = self._current_profiles(position, field, carriers)
        return PNJunctionEquilibrium(
            position=position,
            charge_density=charge_density,
            electrostatic_potential=potential,
            electric_field=field,
            bands=bands,
            carriers=carriers,
            currents=currents,
        )

    @staticmethod
    def _position_array(position: NumpyArrayLike) -> FloatArray:
        positions = np.asarray(position, dtype=np.float64)
        if positions.ndim != 1:
            raise ValueError("position must be one-dimensional")
        if positions.size < 3:
            raise ValueError("position must contain at least three points")
        if not np.all(np.isfinite(positions)):
            raise ValueError("position must contain only finite values")
        if np.any(np.diff(positions) <= 0.0):
            raise ValueError("position must be strictly increasing")
        return positions

    def _solve_poisson(
        self,
        position: FloatArray,
        charge_density: FloatArray,
        *,
        applied_potential: float,
    ) -> tuple[FloatArray, FloatArray]:
        point_count = position.size
        spacing = position[1] - position[0]
        if not np.allclose(np.diff(position), spacing):
            raise ValueError("position must be uniformly spaced")
        matrix = np.zeros((point_count, point_count), dtype=np.float64)
        right_hand_side = np.zeros(point_count, dtype=np.float64)
        for index in range(1, point_count - 1):
            matrix[index, index - 1] = 1.0 / spacing**2
            matrix[index, index] = -2.0 / spacing**2
            matrix[index, index + 1] = 1.0 / spacing**2
            right_hand_side[index] = -charge_density[index] / self.material.permittivity
        matrix[0, 0] = 1.0
        matrix[-1, -1] = 1.0
        right_hand_side[-1] = self.built_in_potential() - applied_potential
        potential = np.linalg.solve(matrix, right_hand_side)
        field = np.zeros_like(potential)
        field[1:-1] = -(potential[2:] - potential[:-2]) / (2.0 * spacing)
        field[0] = field[1]
        field[-1] = field[-2]
        return potential, field

    def _energy_bands(
        self,
        position: FloatArray,
        potential: FloatArray,
    ) -> EnergyBands:
        conduction_raw = -potential
        valence_raw = conduction_raw - self.material.band_gap
        intrinsic_raw = valence_raw + self.material.intrinsic_level_from_valence_band
        intrinsic_target = self.material.thermal_energy * np.log(
            self.doping.acceptor_density / self.material.intrinsic_carrier_density
        )
        shift = intrinsic_target - intrinsic_raw[0]
        return EnergyBands(
            conduction=conduction_raw + shift,
            valence=valence_raw + shift,
            intrinsic=intrinsic_raw + shift,
            fermi=np.zeros_like(position),
        )

    def _carrier_profiles(self, bands: EnergyBands) -> CarrierProfiles:
        intrinsic_density = self.material.intrinsic_carrier_density
        thermal_energy = self.material.thermal_energy
        electron = intrinsic_density * np.exp(
            (bands.fermi - bands.intrinsic) / thermal_energy
        )
        hole = intrinsic_density * np.exp(
            (bands.intrinsic - bands.fermi) / thermal_energy
        )
        return CarrierProfiles(electron=electron, hole=hole)

    def _current_profiles(
        self,
        position: FloatArray,
        field: FloatArray,
        carriers: CarrierProfiles,
    ) -> CurrentProfiles:
        spacing = position[1] - position[0]
        electron_gradient = np.zeros_like(carriers.electron)
        electron_gradient[1:-1] = (
            carriers.electron[2:] - carriers.electron[:-2]
        ) / (2.0 * spacing)
        electron_gradient[0] = electron_gradient[1]
        electron_gradient[-1] = electron_gradient[-2]
        hole_gradient = np.zeros_like(carriers.hole)
        hole_gradient[1:-1] = (
            carriers.hole[2:] - carriers.hole[:-2]
        ) / (2.0 * spacing)
        hole_gradient[0] = hole_gradient[1]
        hole_gradient[-1] = hole_gradient[-2]
        electron_diffusivity = (
            self.material.thermal_energy * self.material.electron_mobility
        )
        hole_diffusivity = self.material.thermal_energy * self.material.hole_mobility
        electron_drift = (
            SI.q * carriers.electron * self.material.electron_mobility * field
        )
        electron_diffusion = SI.q * electron_diffusivity * electron_gradient
        hole_drift = SI.q * carriers.hole * self.material.hole_mobility * field
        hole_diffusion = -SI.q * hole_diffusivity * hole_gradient
        return CurrentProfiles(
            electron_drift=electron_drift,
            electron_diffusion=electron_diffusion,
            electron_total=electron_drift + electron_diffusion,
            hole_drift=hole_drift,
            hole_diffusion=hole_diffusion,
            hole_total=hole_drift + hole_diffusion,
        )
