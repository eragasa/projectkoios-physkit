r"""Unit-aware molar and mass composition for finite species mixtures."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

_MOLE = PhysicalUnit(expression="mole")
_KILOGRAM = PhysicalUnit(expression="kilogram")
_KILOGRAM_PER_MOLE = PhysicalUnit(expression="kilogram / mole")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class MolarMixture(DataObject):
    """Represent ordered species amounts and corresponding molar masses."""

    species: tuple[str, ...]
    mole_amounts: VectorQuantity
    molar_masses: VectorQuantity

    def __post_init__(self) -> None:
        """Check every molar-mixture constructor argument."""
        self._check_arg_species()
        self._check_arg_mole_amounts()
        self._check_arg_molar_masses()

    def _check_arg_species(self) -> None:
        """Require a nonempty ordered tuple of unique nonempty names."""
        if not isinstance(self.species, tuple) or not all(
            isinstance(name, str) for name in self.species
        ):
            raise TypeError("species must be a tuple of strings")
        if not self.species or any(not name for name in self.species):
            raise ValueError("species names must be nonempty")
        if len(set(self.species)) != len(self.species):
            raise ValueError("species names must be unique")

    def _check_arg_mole_amounts(self) -> None:
        """Require finite nonnegative amounts with positive total."""
        if not isinstance(self.mole_amounts, VectorQuantity):
            raise TypeError("mole_amounts must be VectorQuantity")
        if not isinstance(self.mole_amounts.unit, PhysicalUnit):
            raise TypeError("mole_amounts must use a physical amount unit")
        if not ModelSystemUnitConverter.compatible(
            self.mole_amounts.unit,
            _MOLE,
        ):
            raise ValueError("mole_amounts must use units compatible with amount")
        if self.mole_amounts.magnitude.shape != (len(self.species),):
            raise ValueError("mole_amounts must have one value per species")
        if np.any(self.mole_amounts.magnitude < 0.0):
            raise ValueError("mole_amounts must be nonnegative")
        if float(np.sum(self.mole_amounts.magnitude)) <= 0.0:
            raise ValueError("total mole amount must be positive")

    def _check_arg_molar_masses(self) -> None:
        """Require one positive mass-per-amount value per species."""
        if not isinstance(self.molar_masses, VectorQuantity):
            raise TypeError("molar_masses must be VectorQuantity")
        if not isinstance(self.molar_masses.unit, PhysicalUnit):
            raise TypeError("molar_masses must use a physical unit")
        if not ModelSystemUnitConverter.compatible(
            self.molar_masses.unit,
            _KILOGRAM_PER_MOLE,
        ):
            raise ValueError(
                "molar_masses must use units compatible with mass per amount"
            )
        if self.molar_masses.magnitude.shape != (len(self.species),):
            raise ValueError("molar_masses must have one value per species")
        if np.any(self.molar_masses.magnitude <= 0.0):
            raise ValueError("molar_masses must be positive")

    @property
    def total_mole_amount(self) -> ScalarQuantity:
        """Return the total represented amount in the input amount unit."""
        return ScalarQuantity(
            magnitude=float(np.sum(self.mole_amounts.magnitude)),
            unit=self.mole_amounts.unit,
        )

    @property
    def mole_fractions(self) -> VectorQuantity:
        """Return ordered dimensionless mole fractions."""
        return VectorQuantity(
            magnitude=np.asarray(
                self.mole_amounts.magnitude / self.total_mole_amount.magnitude,
                dtype=np.float64,
            ),
            unit=Unitless(),
        )

    @property
    def mean_molar_mass(self) -> ScalarQuantity:
        r"""Return $\overline M=\sum_i x_iM_i$."""
        return ScalarQuantity(
            magnitude=float(
                np.dot(self.mole_fractions.magnitude, self.molar_masses.magnitude)
            ),
            unit=self.molar_masses.unit,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class MassMixture(DataObject):
    """Represent ordered species masses and corresponding molar masses."""

    species: tuple[str, ...]
    mass_amounts: VectorQuantity
    molar_masses: VectorQuantity

    def __post_init__(self) -> None:
        """Check every mass-mixture constructor argument."""
        self._check_arg_species()
        self._check_arg_mass_amounts()
        self._check_arg_molar_masses()

    def _check_arg_species(self) -> None:
        """Require a nonempty ordered tuple of unique nonempty names."""
        if not isinstance(self.species, tuple) or not all(
            isinstance(name, str) for name in self.species
        ):
            raise TypeError("species must be a tuple of strings")
        if not self.species or any(not name for name in self.species):
            raise ValueError("species names must be nonempty")
        if len(set(self.species)) != len(self.species):
            raise ValueError("species names must be unique")

    def _check_arg_mass_amounts(self) -> None:
        """Require finite nonnegative masses with positive total."""
        if not isinstance(self.mass_amounts, VectorQuantity):
            raise TypeError("mass_amounts must be VectorQuantity")
        if not isinstance(self.mass_amounts.unit, PhysicalUnit):
            raise TypeError("mass_amounts must use a physical mass unit")
        if not ModelSystemUnitConverter.compatible(
            self.mass_amounts.unit,
            _KILOGRAM,
        ):
            raise ValueError("mass_amounts must use units compatible with mass")
        if self.mass_amounts.magnitude.shape != (len(self.species),):
            raise ValueError("mass_amounts must have one value per species")
        if np.any(self.mass_amounts.magnitude < 0.0):
            raise ValueError("mass_amounts must be nonnegative")
        if float(np.sum(self.mass_amounts.magnitude)) <= 0.0:
            raise ValueError("total mass amount must be positive")

    def _check_arg_molar_masses(self) -> None:
        """Require one positive mass-per-amount value per species."""
        if not isinstance(self.molar_masses, VectorQuantity):
            raise TypeError("molar_masses must be VectorQuantity")
        if not isinstance(self.molar_masses.unit, PhysicalUnit):
            raise TypeError("molar_masses must use a physical unit")
        if not ModelSystemUnitConverter.compatible(
            self.molar_masses.unit,
            _KILOGRAM_PER_MOLE,
        ):
            raise ValueError(
                "molar_masses must use units compatible with mass per amount"
            )
        if self.molar_masses.magnitude.shape != (len(self.species),):
            raise ValueError("molar_masses must have one value per species")
        if np.any(self.molar_masses.magnitude <= 0.0):
            raise ValueError("molar_masses must be positive")

    @property
    def total_mass_amount(self) -> ScalarQuantity:
        """Return the total represented mass in the input mass unit."""
        return ScalarQuantity(
            magnitude=float(np.sum(self.mass_amounts.magnitude)),
            unit=self.mass_amounts.unit,
        )

    @property
    def mass_fractions(self) -> VectorQuantity:
        """Return ordered dimensionless mass fractions."""
        return VectorQuantity(
            magnitude=np.asarray(
                self.mass_amounts.magnitude / self.total_mass_amount.magnitude,
                dtype=np.float64,
            ),
            unit=Unitless(),
        )

    @property
    def mean_molar_mass(self) -> ScalarQuantity:
        r"""Return $\overline M=(\sum_i w_i/M_i)^{-1}$."""
        reciprocal = float(
            np.sum(self.mass_fractions.magnitude / self.molar_masses.magnitude)
        )
        return ScalarQuantity(
            magnitude=1.0 / reciprocal,
            unit=self.molar_masses.unit,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class MolarToMassMixtureConversionRequest(DataObject):
    """Request conversion from molar amounts to masses."""

    mixture: MolarMixture
    target_mass_unit: PhysicalUnit

    def __post_init__(self) -> None:
        """Check every molar-to-mass request argument."""
        self._check_arg_mixture()
        self._check_arg_target_mass_unit()

    def _check_arg_mixture(self) -> None:
        """Require a molar mixture."""
        if not isinstance(self.mixture, MolarMixture):
            raise TypeError("mixture must be MolarMixture")

    def _check_arg_target_mass_unit(self) -> None:
        """Require a physical mass unit."""
        if not isinstance(self.target_mass_unit, PhysicalUnit):
            raise TypeError("target_mass_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.target_mass_unit,
            _KILOGRAM,
        ):
            raise ValueError("target_mass_unit must be compatible with mass")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class MolarToMassMixtureConversion(ResultsObject):
    """Correlate a molar-to-mass request with its mass mixture."""

    request: MolarToMassMixtureConversionRequest
    mixture: MassMixture

    def __post_init__(self) -> None:
        """Check every molar-to-mass response argument."""
        self._check_arg_request()
        self._check_arg_mixture()

    def _check_arg_request(self) -> None:
        """Require the complete conversion request."""
        if not isinstance(self.request, MolarToMassMixtureConversionRequest):
            raise TypeError("request must be MolarToMassMixtureConversionRequest")

    def _check_arg_mixture(self) -> None:
        """Require a mass mixture with matching species order."""
        if not isinstance(self.mixture, MassMixture):
            raise TypeError("mixture must be MassMixture")
        if self.mixture.species != self.request.mixture.species:
            raise ValueError("converted mixture must preserve species order")


@dataclass(frozen=True, slots=True, kw_only=True)
class MassToMolarMixtureConversionRequest(DataObject):
    """Request conversion from masses to molar amounts."""

    mixture: MassMixture
    target_amount_unit: PhysicalUnit

    def __post_init__(self) -> None:
        """Check every mass-to-molar request argument."""
        self._check_arg_mixture()
        self._check_arg_target_amount_unit()

    def _check_arg_mixture(self) -> None:
        """Require a mass mixture."""
        if not isinstance(self.mixture, MassMixture):
            raise TypeError("mixture must be MassMixture")

    def _check_arg_target_amount_unit(self) -> None:
        """Require a physical amount unit."""
        if not isinstance(self.target_amount_unit, PhysicalUnit):
            raise TypeError("target_amount_unit must be PhysicalUnit")
        if not ModelSystemUnitConverter.compatible(
            self.target_amount_unit,
            _MOLE,
        ):
            raise ValueError("target_amount_unit must be compatible with amount")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class MassToMolarMixtureConversion(ResultsObject):
    """Correlate a mass-to-molar request with its molar mixture."""

    request: MassToMolarMixtureConversionRequest
    mixture: MolarMixture

    def __post_init__(self) -> None:
        """Check every mass-to-molar response argument."""
        self._check_arg_request()
        self._check_arg_mixture()

    def _check_arg_request(self) -> None:
        """Require the complete conversion request."""
        if not isinstance(self.request, MassToMolarMixtureConversionRequest):
            raise TypeError("request must be MassToMolarMixtureConversionRequest")

    def _check_arg_mixture(self) -> None:
        """Require a molar mixture with matching species order."""
        if not isinstance(self.mixture, MolarMixture):
            raise TypeError("mixture must be MolarMixture")
        if self.mixture.species != self.request.mixture.species:
            raise ValueError("converted mixture must preserve species order")


@dataclass(frozen=True, slots=True)
class MixtureCompositionConverter:
    """Convert ordered finite mixtures between molar and mass amounts."""

    def convert_to_mass(
        self,
        *,
        request: MolarToMassMixtureConversionRequest,
    ) -> MolarToMassMixtureConversion:
        """Convert species mole amounts to masses in the requested unit."""
        if not isinstance(request, MolarToMassMixtureConversionRequest):
            raise TypeError("request must be MolarToMassMixtureConversionRequest")
        mole_amounts = ModelSystemUnitConverter.convert_vector(
            request.mixture.mole_amounts,
            _MOLE,
        )
        molar_masses = ModelSystemUnitConverter.convert_vector(
            request.mixture.molar_masses,
            _KILOGRAM_PER_MOLE,
        )
        masses = ModelSystemUnitConverter.convert_vector(
            VectorQuantity(
                magnitude=np.asarray(
                    mole_amounts.magnitude * molar_masses.magnitude,
                    dtype=np.float64,
                ),
                unit=_KILOGRAM,
            ),
            request.target_mass_unit,
        )
        return MolarToMassMixtureConversion(
            request=request,
            mixture=MassMixture(
                species=request.mixture.species,
                mass_amounts=masses,
                molar_masses=request.mixture.molar_masses,
            ),
        )

    def convert_to_molar(
        self,
        *,
        request: MassToMolarMixtureConversionRequest,
    ) -> MassToMolarMixtureConversion:
        """Convert species masses to mole amounts in the requested unit."""
        if not isinstance(request, MassToMolarMixtureConversionRequest):
            raise TypeError("request must be MassToMolarMixtureConversionRequest")
        masses = ModelSystemUnitConverter.convert_vector(
            request.mixture.mass_amounts,
            _KILOGRAM,
        )
        molar_masses = ModelSystemUnitConverter.convert_vector(
            request.mixture.molar_masses,
            _KILOGRAM_PER_MOLE,
        )
        mole_amounts = ModelSystemUnitConverter.convert_vector(
            VectorQuantity(
                magnitude=np.asarray(
                    masses.magnitude / molar_masses.magnitude,
                    dtype=np.float64,
                ),
                unit=_MOLE,
            ),
            request.target_amount_unit,
        )
        return MassToMolarMixtureConversion(
            request=request,
            mixture=MolarMixture(
                species=request.mixture.species,
                mole_amounts=mole_amounts,
                molar_masses=request.mixture.molar_masses,
            ),
        )
