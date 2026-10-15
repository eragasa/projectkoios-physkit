# physkit/units/systems.py
#
# Author: Eugene Joseph M. Ragasa
# Affiliation: De La Salle University
#
# Canonical unit-system views for physkit quantities.
#
# These classes provide *conventional defaults* only.
# They do NOT perform conversions and carry no state.

from enum import StrEnum

from .charge import Charge
from .density import Density
from .dipole import Dipole
from .electricfield import ElectricField
from .energy import Energy
from .force import Force
from .length import Length
from .mass import Mass
from .pressure import Pressure
from .quantities import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    ModelSystemUnit,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
)
from .temperature import Temperature
from .time import Time
from .torque import Torque
from .velocity import Velocity
from .viscosity import Viscosity


class UnitSystem(StrEnum):
    """Select a coherent numerical scale for model-system calculations."""

    NONDIMENSIONAL = "nondimensional"
    SI = "si"
    METAL = "metal"

    @property
    def length_unit(self) -> ModelSystemUnit:
        """Return the system's numerical length unit."""
        if self is UnitSystem.NONDIMENSIONAL:
            return Unitless()
        if self is UnitSystem.SI:
            return PhysicalUnit("meter")
        return PhysicalUnit("angstrom")

    @property
    def mass_unit(self) -> ModelSystemUnit:
        """Return the system's numerical particle-mass unit.

        Metal-style mass values use daltons. Their magnitudes match the
        conventional grams-per-mole values while retaining particle-mass
        dimensions suitable for quantum mechanics.
        """
        if self is UnitSystem.NONDIMENSIONAL:
            return Unitless()
        if self is UnitSystem.SI:
            return PhysicalUnit("kilogram")
        return PhysicalUnit("dalton")

    @property
    def time_unit(self) -> ModelSystemUnit:
        """Return the system's numerical time unit."""
        if self is UnitSystem.NONDIMENSIONAL:
            return Unitless()
        if self is UnitSystem.SI:
            return PhysicalUnit("second")
        return PhysicalUnit("picosecond")

    @property
    def energy_unit(self) -> ModelSystemUnit:
        """Return the system's numerical energy unit."""
        if self is UnitSystem.NONDIMENSIONAL:
            return Unitless()
        if self is UnitSystem.SI:
            return PhysicalUnit("joule")
        return PhysicalUnit("electron_volt")

    @property
    def action_unit(self) -> ModelSystemUnit:
        """Return the system's numerical action unit."""
        if self is UnitSystem.NONDIMENSIONAL:
            return Unitless()
        if self is UnitSystem.SI:
            return PhysicalUnit("joule * second")
        return PhysicalUnit("electron_volt * picosecond")

    @property
    def hbar(self) -> ScalarQuantity:
        """Return reduced Planck's constant in the selected action unit."""
        if self is UnitSystem.NONDIMENSIONAL:
            return ScalarQuantity(1.0, Unitless())
        hbar_si = ScalarQuantity(
            1.054_571_817e-34,
            PhysicalUnit("joule * second"),
        )
        return MODEL_SYSTEM_UNIT_CONVERTER.convert_scalar(
            hbar_si,
            self.action_unit,
        )


class UnitsSI:
    """
    SI unit system (absolute).

    Canonical mechanical + EM SI units.
    """

    name = "SI"

    # --- kinematics ---
    length = Length.Units.m
    time = Time.Units.s
    velocity = Velocity.Units.m_per_s

    # --- mechanics ---
    mass = Mass.Units.kg
    force = Force.Units.N
    energy = Energy.Units.J
    torque = Torque.Units.N_m
    pressure = Pressure.Units.Pa
    viscosity = Viscosity.Units.Pa_s  # dynamic viscosity μ

    # --- thermodynamics ---
    temperature = Temperature.Units.K
    density = Density.Units.kg_per_m_dim  # m^dim (2D/3D aware)

    # --- electromagnetism ---
    charge = Charge.Units.C
    dipole = Dipole.Units.C_m
    electric_field = ElectricField.Units.V_per_m


class UnitsCGS:
    """
    CGS unit system (cm–g–s, electrostatic/esu).

    Mechanical CGS + electrostatic CGS (Gaussian-compatible).
    """

    name = "CGS"

    # --- kinematics ---
    length = Length.Units.cm
    time = Time.Units.s
    velocity = Velocity.Units.cm_per_s

    # --- mechanics ---
    mass = Mass.Units.g
    force = Force.Units.dyn
    energy = Energy.Units.erg
    torque = Torque.Units.dyn_cm
    pressure = Pressure.Units.Ba  # barye = dyn/cm^2
    viscosity = Viscosity.Units.Poise  # dynamic viscosity μ

    # --- thermodynamics ---
    temperature = Temperature.Units.K
    density = Density.Units.g_per_cm_dim  # cm^dim aware

    # --- electromagnetism (CGS–esu) ---
    charge = Charge.Units.esu  # statcoulomb
    dipole = Dipole.Units.esu_cm
    electric_field = ElectricField.Units.statV_per_cm
    # equivalent to dyn / esu


class UnitsElectron:
    """
    Electron-scale unit system (hybrid atomic units).

    Electronic energies and lengths are atomic-scale,
    while time, temperature, pressure, and fields remain
    in laboratory units.
    """

    name = "lammps.electron"

    # --- kinematics ---
    length = Length.Units.bohr
    time = Time.Units.fs
    velocity = Velocity.Units.bohr_per_t0  # atomic velocity unit

    # --- mechanics ---
    mass = Mass.Units.amu
    energy = Energy.Units.Ha
    force = Force.Units.Ha_per_bohr

    # --- thermodynamics ---
    temperature = Temperature.Units.K
    pressure = Pressure.Units.Pa

    # --- electromagnetism ---
    charge = Charge.Units.e  # dimensionless multiples of e
    dipole = Dipole.Units.Debye
    electric_field = ElectricField.Units.V_per_cm


class UnitsHartree:
    """
    Hartree atomic units (a.u.).

    Fundamental constants set to unity:
        ħ = 1
        m_e = 1
        e = 1
        4πϵ₀ = 1

    This is a Hamiltonian-normalized unit system,
    not a laboratory unit system.
    """

    name = "hartree"

    # --- mechanics / quantum ---
    length = Length.Units.bohr  # a₀
    mass = Mass.Units.me  # electron mass
    time = Time.Units.atomic  # t₀ = ħ / E_h
    energy = Energy.Units.Ha  # Hartree
    velocity = Velocity.Units.atomic  # a₀ / t₀
    force = Force.Units.Ha_per_bohr

    # --- electromagnetism ---
    charge = Charge.Units.atomic  # dimensionless (±1)
    electric_field = ElectricField.Units.atomic

    # --- derived (conceptual, not SI) ---
    pressure = Pressure.Units.atomic  # E_h / a₀³
    density = Density.Units.atomic  # m_e / a₀^dim


class UnitsReal:
    """
    Real unit system (LAMMPS 'real')

    https://docs.lammps.org/units.html
    """

    name = "lammps.real"

    mass = Mass.Units.g_per_mol
    length = Length.Units.A
    time = Time.Units.fs
    energy = Energy.Units.kcal_per_mol

    velocity = Velocity.Units.A_per_fs
    force = Force.Units.kcal_per_mol_A
    torque = Torque.Units.kcal_per_mol

    temperature = Temperature.Units.K
    pressure = Pressure.Units.atm

    viscosity = Viscosity.Units.Poise  # dynamic viscosity

    charge = Charge.Units.e
    dipole = Dipole.Units.e_A

    electric_field = ElectricField.Units.V_per_A

    density = Density.Units.g_per_cm3  # (3D); in general LAMMPS uses g/cm^dim


class UnitsMetal:
    """
    Metal / condensed-matter unit system.

    Practical MD units for crystalline solids and metals.
    """

    name = "lammps.metal"

    # --- kinematics ---
    length = Length.Units.A
    time = Time.Units.ps
    velocity = Velocity.Units.A_per_ps

    # --- mechanics ---
    mass = Mass.Units.g_per_mol
    energy = Energy.Units.eV
    force = Force.Units.eV_per_A
    torque = Torque.Units.eV
    pressure = Pressure.Units.bar
    viscosity = Viscosity.Units.Poise  # dynamic viscosity μ

    # --- thermodynamics ---
    temperature = Temperature.Units.K
    density = Density.Units.g_per_cm_dim  # cm^dim aware

    # --- electromagnetism ---
    charge = Charge.Units.e  # ±1 = ±e
    dipole = Dipole.Units.e_A
    electric_field = ElectricField.Units.V_per_A


class UnitsImperial:
    """
    Imperial unit system (absolute).

    Uses absolute force (lbf), not mass-based units.
    """

    length = Length.Units.ft
    mass = Mass.Units.lbm
    force = Force.Units.lbf
    pressure = Pressure.Units.psi
    temperature = Temperature.Units.R


class UnitsUSCS:
    """
    US Customary System (engineering, absolute)
    """

    length = Length.Units.ft
    force = Force.Units.lbf
    pressure = Pressure.Units.psi
    temperature = Temperature.Units.R
