from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from projectkoios.physkit.constants import ConstantsSI


class MolarMassUnit(str, Enum):  # noqa: UP042 - Preserve legacy Enum string behavior.
    G_PER_MOL = "g/mol"
    KG_PER_MOL = "kg/mol"
    KG_PER_KMOL = "kg/kmol"  # common in chem eng


@dataclass(frozen=True)
class MolarMass:
    """
    Canonical storage: kg/mol.
    """

    kg_per_mol: float

    def __post_init__(self) -> None:
        if self.kg_per_mol <= 0.0:
            raise ValueError("MolarMass must be > 0.")

    @staticmethod
    def from_value(value: float, unit: MolarMassUnit) -> MolarMass:
        if value <= 0.0:
            raise ValueError("value must be > 0.")
        if unit == MolarMassUnit.KG_PER_MOL:
            kg_per_mol = value
        elif unit == MolarMassUnit.G_PER_MOL:
            kg_per_mol = value * 1e-3
        elif unit == MolarMassUnit.KG_PER_KMOL:
            kg_per_mol = value * 1e-3  # 1 kg/kmol = 1e-3 kg/mol
        else:
            raise ValueError(f"Unsupported unit: {unit}")
        return MolarMass(kg_per_mol)

    def to(self, unit: MolarMassUnit) -> float:
        if unit == MolarMassUnit.KG_PER_MOL:
            return self.kg_per_mol
        elif unit == MolarMassUnit.G_PER_MOL or unit == MolarMassUnit.KG_PER_KMOL:
            return self.kg_per_mol * 1e3
        else:
            raise ValueError(f"Unsupported unit: {unit}")

    def to_particle_mass(self) -> ParticleMass:
        return ParticleMass(self.kg_per_mol / ConstantsSI.N_A)


class ParticleMassUnit(str, Enum):  # noqa: UP042 - Preserve legacy Enum string behavior.
    KG = "kg"
    AMU = "amu"  # optional convenience


@dataclass(frozen=True)
class ParticleMass:
    """
    Canonical storage: kg (per particle).
    """

    kg: float

    def __post_init__(self) -> None:
        if self.kg <= 0.0:
            raise ValueError("ParticleMass must be > 0.")

    @staticmethod
    def from_value(value: float, unit: ParticleMassUnit) -> ParticleMass:
        if value <= 0.0:
            raise ValueError("value must be > 0.")
        if unit == ParticleMassUnit.KG:
            kg = value
        elif unit == ParticleMassUnit.AMU:
            kg = value * ConstantsSI.atomic_mass_constant  # if you have it
        else:
            raise ValueError(f"Unsupported unit: {unit}")
        return ParticleMass(kg)

    def to(self, unit: ParticleMassUnit) -> float:
        if unit == ParticleMassUnit.KG:
            return self.kg
        elif unit == ParticleMassUnit.AMU:
            return self.kg / ConstantsSI.atomic_mass_constant
        else:
            raise ValueError(f"Unsupported unit: {unit}")

    def to_molar_mass(self) -> MolarMass:
        return MolarMass(self.kg * ConstantsSI.N_A)
