# src/python/projectkoios/physkit/constants.py
# Eugene Joseph M. Ragasa 2025-2026

r"""
Physical constants for explicit numerical unit systems.

The module provides immutable constants containers for direct numerical
work. A physical model receives a constants object explicitly, and all
model parameters must be expressed consistently in that constants
object's unit system.

Examples
--------
Construct a model using SI constants:

>>> model = Piab1D(
...     x_lower=0.0,
...     x_upper=1.0e-9,
...     mass=SI.me0,
...     constants=SI,
... )

Construct the equivalent model using Gaussian CGS constants:

>>> model = Piab1D(
...     x_lower=0.0,
...     x_upper=1.0e-7,
...     mass=GAUSSIAN_CGS.me0,
...     constants=GAUSSIAN_CGS,
... )

Unit consistency is the responsibility of the model configuration. This
module does not perform dimensional analysis, automatic unit conversion,
or conversion between SI and Gaussian CGS equation forms.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import (
    Final,
    Protocol,
    runtime_checkable,
)


@runtime_checkable
class PhysicalConstantsProtocol(Protocol):
    """Structural interface shared by numerical constants containers."""

    a0: float
    q: float
    k_B: float
    eps0: float | None
    me0: float
    N_A: float
    R_g: float
    h: float
    hbar: float
    m_u: float
    m_u_u: float | None
    c: float


@dataclass(
    frozen=True,
    slots=True,
)
class ConstantsSI:
    """
    Physical constants in SI units.

    Attributes
    ----------
    a0:
        Bohr radius in meters.
    q:
        Elementary charge in coulombs.
    k_B:
        Boltzmann constant in joules per kelvin.
    eps0:
        Vacuum permittivity in farads per meter.
    me0:
        Electron rest mass in kilograms.
    N_A:
        Avogadro constant in reciprocal moles.
    R_g:
        Molar gas constant in joules per mole kelvin.
    h:
        Planck constant in joule seconds.
    hbar:
        Reduced Planck constant in joule seconds.
    m_u:
        Atomic mass constant in kilograms.
    m_u_u:
        Standard uncertainty of the atomic mass constant in kilograms.
    c:
        Speed of light in meters per second.
    """

    a0: float = 5.291_772_109_03e-11
    q: float = 1.602_176_634e-19
    k_B: float = 1.380_649e-23
    eps0: float = 8.854_187_812_8e-12
    me0: float = 9.109_383_713_9e-31
    N_A: float = 6.022_140_76e23
    R_g: float = 8.314_462_618_153_24
    h: float = 6.626_070_15e-34
    hbar: float = h / (2.0 * math.pi)
    m_u: float = 1.660_539_068_92e-27
    m_u_u: float = 0.000_000_000_52e-27
    c: float = 299_792_458.0


@dataclass(
    frozen=True,
    slots=True,
)
class ConstantsGaussianCGS:
    """
    Physical constants in Gaussian CGS units.

    Attributes
    ----------
    a0:
        Bohr radius in centimeters.
    q:
        Elementary charge in statcoulombs.
    k_B:
        Boltzmann constant in ergs per kelvin.
    me0:
        Electron rest mass in grams.
    N_A:
        Avogadro constant in reciprocal moles.
    R_g:
        Molar gas constant in ergs per mole kelvin.
    h:
        Planck constant in erg seconds.
    hbar:
        Reduced Planck constant in erg seconds.
    m_u:
        Atomic mass constant in grams.
    m_u_u:
        Standard uncertainty of the atomic mass constant in grams.
    c:
        Speed of light in centimeters per second.

    Notes
    -----
    Gaussian CGS electromagnetic equations do not use vacuum permittivity in
    the same form as SI electromagnetic equations. ``eps0`` is therefore
    represented explicitly as ``None`` rather than as a fabricated numerical
    value.
    """

    a0: float = 5.291_772_109_03e-9
    q: float = 4.803_204_712_57e-10
    k_B: float = 1.380_649e-16
    eps0: None = None
    me0: float = 9.109_383_713_9e-28
    N_A: float = 6.022_140_76e23
    R_g: float = 8.314_462_618_153_24e7
    h: float = 6.626_070_15e-27
    hbar: float = h / (2.0 * math.pi)
    m_u: float = 1.660_539_068_92e-24
    m_u_u: float = 0.000_000_000_52e-24
    c: float = 2.997_924_58e10


ConstantsCGS = ConstantsGaussianCGS

type Constants = ConstantsSI | ConstantsGaussianCGS


SI: Final[ConstantsSI] = ConstantsSI()

GAUSSIAN_CGS: Final[ConstantsGaussianCGS] = ConstantsGaussianCGS()
