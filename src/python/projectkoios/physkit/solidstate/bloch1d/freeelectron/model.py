r"""Exact folded free-electron Bloch spectrum in one dimension."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.constants import SI
from projectkoios.physkit.core.actions import ActionizedDataObject, DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.solidstate.bloch1d.cell import (
    BlochPrimitiveCell1D,
)
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    PhysicalUnit,
    ScalarQuantity,
    VectorQuantity,
)

ModeIndex1D = tuple[int]


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class FreeElectronBloch1DSpectrumRequest(DataObject):
    """Request folded free-electron modes at one Bloch wave number."""

    model: FreeElectronBloch1D
    mode_indices: tuple[ModeIndex1D, ...]
    bloch_wave_number: ScalarQuantity

    def __post_init__(self) -> None:
        """Check every spectrum-request argument."""
        self._check_arg_model()
        self._check_arg_mode_indices()
        self._check_arg_bloch_wave_number()

    def _check_arg_model(self) -> None:
        if not isinstance(self.model, FreeElectronBloch1D):
            raise TypeError("model must be FreeElectronBloch1D")

    def _check_arg_mode_indices(self) -> None:
        if type(self.mode_indices) is not tuple:
            raise TypeError("mode_indices must be a tuple")
        if len(self.mode_indices) == 0:
            raise ValueError("mode_indices must not be empty")
        for mode in self.mode_indices:
            if type(mode) is not tuple or len(mode) != 1:
                raise TypeError("each mode index must be a one-integer tuple")
            if type(mode[0]) is not int:
                raise TypeError("mode-index components must be built-in integers")
        if len(set(self.mode_indices)) != len(self.mode_indices):
            raise ValueError("mode_indices must be unique")

    def _check_arg_bloch_wave_number(self) -> None:
        if not isinstance(self.bloch_wave_number, ScalarQuantity):
            raise TypeError("bloch_wave_number must be ScalarQuantity")
        if self.bloch_wave_number.unit != self.model.cell.reciprocal_basis.unit:
            raise ValueError(
                "bloch_wave_number must use the inverse selected length unit"
            )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class FreeElectronBloch1DSpectrum(ResultsObject):
    """Return physical wave numbers and energies in request order."""

    request: FreeElectronBloch1DSpectrumRequest
    wave_numbers: VectorQuantity
    energies: VectorQuantity

    def __post_init__(self) -> None:
        """Check every spectrum-result argument."""
        self._check_arg_request()
        self._check_arg_wave_numbers()
        self._check_arg_energies()

    def _check_arg_request(self) -> None:
        if not isinstance(self.request, FreeElectronBloch1DSpectrumRequest):
            raise TypeError("request must be FreeElectronBloch1DSpectrumRequest")

    def _check_arg_wave_numbers(self) -> None:
        if not isinstance(self.wave_numbers, VectorQuantity):
            raise TypeError("wave_numbers must be VectorQuantity")
        if self.wave_numbers.unit != self.request.model.cell.reciprocal_basis.unit:
            raise ValueError("wave_numbers must use the reciprocal-length unit")
        if self.wave_numbers.magnitude.shape != (len(self.request.mode_indices),):
            raise ValueError("wave_numbers shape must match requested modes")

    def _check_arg_energies(self) -> None:
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if self.energies.unit != self.request.model.cell.unit_system.energy_unit:
            raise ValueError("energies must use the selected energy unit")
        if self.energies.magnitude.shape != (len(self.request.mode_indices),):
            raise ValueError("energies shape must match requested modes")
        if np.any(self.energies.magnitude < 0.0):
            raise ValueError("energies must be nonnegative")


class FreeElectronBloch1DSpectrumEvaluator(
    DataObjectActionizer[
        FreeElectronBloch1DSpectrumRequest,
        FreeElectronBloch1DSpectrum,
    ]
):
    """Evaluate exact folded bands using the bare electron mass."""

    def action(
        self,
        *,
        request: FreeElectronBloch1DSpectrumRequest,
    ) -> FreeElectronBloch1DSpectrum:
        r"""Return $q_n=k+2\pi n/a$ and its exact free-electron energy."""
        if not isinstance(request, FreeElectronBloch1DSpectrumRequest):
            raise TypeError("request must be FreeElectronBloch1DSpectrumRequest")
        model = request.model
        reciprocal_basis = model.cell.reciprocal_basis
        indices = np.asarray(
            [mode[0] for mode in request.mode_indices],
            dtype=np.float64,
        )
        wave_numbers = (
            request.bloch_wave_number.magnitude + reciprocal_basis.magnitude * indices
        )
        mass = model.electron_mass
        hbar = model.cell.unit_system.hbar
        represented_energy_unit = PhysicalUnit(
            f"({hbar.unit.expression}) ** 2 / "
            f"(({mass.unit.expression}) * "
            f"({model.cell.length.unit.expression}) ** 2)"
        )
        conversion_factor = MODEL_SYSTEM_UNIT_CONVERTER.conversion_factor(
            represented_energy_unit,
            model.cell.unit_system.energy_unit,
        )
        energies = (
            0.5
            * hbar.magnitude**2
            * wave_numbers**2
            / mass.magnitude
            * conversion_factor
        )
        return FreeElectronBloch1DSpectrum(
            request=request,
            wave_numbers=VectorQuantity(
                magnitude=np.asarray(wave_numbers, dtype=np.float64),
                unit=reciprocal_basis.unit,
            ),
            energies=VectorQuantity(
                magnitude=np.asarray(energies, dtype=np.float64),
                unit=model.cell.unit_system.energy_unit,
            ),
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class FreeElectronBloch1D(
    ActionizedDataObject[
        FreeElectronBloch1DSpectrumRequest,
        FreeElectronBloch1DSpectrum,
    ]
):
    """Define bare free-electron Bloch fibers on one primitive cell."""

    cell: BlochPrimitiveCell1D
    actionizer: DataObjectActionizer[
        FreeElectronBloch1DSpectrumRequest,
        FreeElectronBloch1DSpectrum,
    ] = field(
        default_factory=FreeElectronBloch1DSpectrumEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the primitive-cell argument."""
        self._check_arg_cell()

    def _check_arg_cell(self) -> None:
        if not isinstance(self.cell, BlochPrimitiveCell1D):
            raise TypeError("cell must be BlochPrimitiveCell1D")

    @property
    def electron_mass(self) -> ScalarQuantity:
        """Return the bare electron mass in the selected unit system."""
        return MODEL_SYSTEM_UNIT_CONVERTER.convert_scalar(
            ScalarQuantity(SI.me0, PhysicalUnit("kilogram")),
            self.cell.unit_system.mass_unit,
        )

    def evaluate_modes(
        self,
        *,
        mode_indices: tuple[ModeIndex1D, ...],
        bloch_wave_number: ScalarQuantity,
    ) -> FreeElectronBloch1DSpectrum:
        """Evaluate folded modes at one Bloch wave number."""
        return self._respond(
            request=FreeElectronBloch1DSpectrumRequest(
                model=self,
                mode_indices=mode_indices,
                bloch_wave_number=bloch_wave_number,
            )
        )
