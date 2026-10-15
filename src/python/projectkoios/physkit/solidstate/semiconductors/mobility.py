"""Empirical carrier-mobility models with explicit conventional units."""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np

from projectkoios.numerics.numpy.typing import Float64Array, NumpyArrayLike
from projectkoios.physkit.constants import SI


def _positive_array(value: NumpyArrayLike, *, name: str) -> Float64Array:
    array = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    if np.any(array <= 0.0):
        raise ValueError(f"{name} must be positive")
    return array


def _positive_scalar(value: float, *, name: str) -> None:
    if type(value) not in (int, float):
        raise TypeError(f"{name} must be a built-in real number")
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")


@dataclass(frozen=True, slots=True, kw_only=True)
class CarrierMobilityModel:
    """Represent lattice and ionized-impurity mobility for one carrier.

    Effective mass is measured in kilograms, temperature in kelvin, impurity
    concentration in inverse cubic centimeters, and mobility in square
    centimeters per volt second.
    """

    effective_mass: float
    lattice_mobility_reference: float
    impurity_mobility_reference: float
    impurity_concentration_reference: float
    lattice_temperature_reference: float = 300.0
    impurity_temperature_reference: float = 300.0
    lattice_exponent: float = 1.5
    impurity_exponent: float = 1.5

    def __post_init__(self) -> None:
        """Check every empirical-model parameter."""
        for name in (
            "effective_mass",
            "lattice_mobility_reference",
            "impurity_mobility_reference",
            "impurity_concentration_reference",
            "lattice_temperature_reference",
            "impurity_temperature_reference",
            "lattice_exponent",
            "impurity_exponent",
        ):
            _positive_scalar(getattr(self, name), name=name)

    def thermal_velocity(self, temperature: NumpyArrayLike) -> Float64Array:
        """Return the three-dimensional root-mean-square thermal velocity."""
        temperature_array = _positive_array(temperature, name="temperature")
        return np.sqrt(3.0 * SI.k_B * temperature_array / self.effective_mass)

    def relaxation_time(self, mobility: NumpyArrayLike) -> Float64Array:
        """Convert mobility in cm²/(V s) to relaxation time in seconds."""
        mobility_array = _positive_array(mobility, name="mobility")
        mobility_si = mobility_array * 1.0e-4
        return mobility_si * self.effective_mass / SI.q

    def mobility_from_relaxation_time(
        self,
        relaxation_time: NumpyArrayLike,
    ) -> Float64Array:
        """Convert relaxation time in seconds to mobility in cm²/(V s)."""
        time_array = _positive_array(relaxation_time, name="relaxation_time")
        mobility_si = SI.q * time_array / self.effective_mass
        return mobility_si * 1.0e4

    def lattice_mobility(self, temperature: NumpyArrayLike) -> Float64Array:
        """Return the acoustic-phonon power-law mobility contribution."""
        temperature_array = _positive_array(temperature, name="temperature")
        return self.lattice_mobility_reference * (
            temperature_array / self.lattice_temperature_reference
        ) ** (-self.lattice_exponent)

    def impurity_mobility(
        self,
        temperature: NumpyArrayLike,
        impurity_concentration: NumpyArrayLike,
    ) -> Float64Array:
        """Return the screened ionized-impurity mobility contribution."""
        temperature_array = _positive_array(temperature, name="temperature")
        concentration_array = _positive_array(
            impurity_concentration,
            name="impurity_concentration",
        )
        return (
            self.impurity_mobility_reference
            * (temperature_array / self.impurity_temperature_reference)
            ** self.impurity_exponent
            * (self.impurity_concentration_reference / concentration_array)
        )

    def total_mobility(
        self,
        temperature: NumpyArrayLike,
        impurity_concentration: NumpyArrayLike,
    ) -> Float64Array:
        """Combine scattering channels using Matthiessen's rule."""
        lattice = self.lattice_mobility(temperature)
        impurity = self.impurity_mobility(
            temperature,
            impurity_concentration,
        )
        return 1.0 / (1.0 / lattice + 1.0 / impurity)

    def drift_velocity(
        self,
        temperature: NumpyArrayLike,
        impurity_concentration: NumpyArrayLike,
        electric_field: NumpyArrayLike,
    ) -> Float64Array:
        """Return mobility magnitude times electric field in meters per second."""
        field_array = np.asarray(electric_field, dtype=np.float64)
        if not np.all(np.isfinite(field_array)):
            raise ValueError("electric_field must contain only finite values")
        mobility_si = self.total_mobility(temperature, impurity_concentration) * 1.0e-4
        return mobility_si * field_array


@dataclass(frozen=True, slots=True, kw_only=True)
class SemiconductorMobilityModel:
    """Own electron and hole empirical mobility models."""

    electron: CarrierMobilityModel
    hole: CarrierMobilityModel

    def __post_init__(self) -> None:
        """Check both carrier-model arguments."""
        if not isinstance(self.electron, CarrierMobilityModel):
            raise TypeError("electron must be CarrierMobilityModel")
        if not isinstance(self.hole, CarrierMobilityModel):
            raise TypeError("hole must be CarrierMobilityModel")


@dataclass(frozen=True, slots=True, kw_only=True)
class CarrierTrajectory:
    """Store one deterministic realization of a one-dimensional trajectory."""

    time: Float64Array
    position: Float64Array
    velocity: Float64Array


def simulate_electron_trajectory(
    *,
    model: CarrierMobilityModel,
    temperature: float = 300.0,
    impurity_concentration: float = 1.0e16,
    electric_field: float = 1.0e3,
    time_step: float = 1.0e-14,
    step_count: int = 5_000,
    seed: int | None = 1_234,
) -> CarrierTrajectory:
    """Simulate the notebook's relaxation-time electron random walk.

    The model uses a first-order collision probability ``time_step / tau``.
    It is a pedagogical stochastic trajectory, not a transport acceptance
    model or a self-consistent Boltzmann solver.
    """
    if not isinstance(model, CarrierMobilityModel):
        raise TypeError("model must be CarrierMobilityModel")
    _positive_scalar(temperature, name="temperature")
    _positive_scalar(impurity_concentration, name="impurity_concentration")
    _positive_scalar(time_step, name="time_step")
    if type(electric_field) not in (int, float):
        raise TypeError("electric_field must be a built-in real number")
    if not math.isfinite(electric_field):
        raise ValueError("electric_field must be finite")
    if type(step_count) is not int:
        raise TypeError("step_count must be a built-in integer")
    if step_count <= 0:
        raise ValueError("step_count must be positive")
    if seed is not None and type(seed) is not int:
        raise TypeError("seed must be a built-in integer or None")

    mobility = float(model.total_mobility(temperature, impurity_concentration))
    relaxation_time = float(model.relaxation_time(mobility))
    collision_probability = time_step / relaxation_time
    if collision_probability > 1.0:
        raise ValueError("time_step must not exceed the relaxation time")
    if collision_probability > 0.1:
        warnings.warn(
            "time_step is not much smaller than the relaxation time",
            RuntimeWarning,
            stacklevel=2,
        )

    thermal_velocity = float(model.thermal_velocity(temperature))
    component_standard_deviation = thermal_velocity / np.sqrt(3.0)
    acceleration = -SI.q * electric_field / model.effective_mass
    generator = np.random.default_rng(seed)

    time = np.linspace(0.0, step_count * time_step, step_count + 1)
    position = np.zeros(step_count + 1, dtype=np.float64)
    velocity = np.zeros(step_count + 1, dtype=np.float64)
    velocity[0] = generator.normal(0.0, component_standard_deviation)

    for step in range(step_count):
        if generator.random() < collision_probability:
            velocity[step] = generator.normal(0.0, component_standard_deviation)
        velocity[step + 1] = velocity[step] + acceleration * time_step
        position[step + 1] = position[step] + velocity[step + 1] * time_step

    return CarrierTrajectory(
        time=time,
        position=position,
        velocity=velocity,
    )
