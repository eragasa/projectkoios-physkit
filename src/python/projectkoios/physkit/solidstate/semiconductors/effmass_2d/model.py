r"""Bloch-mode spectra for a two-dimensional periodic primitive cell."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject, DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.cell import (
    PeriodicPrimitiveCell2D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.tensor import (
    CartesianEffectiveMassTensor2D,
)
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    MatrixQuantity,
    PhysicalUnit,
    VectorQuantity,
)

ModeIndex2D = tuple[int, int]


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicFreeParticle2DSpectrumRequest(DataObject):
    """Request exact kinetic energies for identified reciprocal modes."""

    model: PeriodicFreeParticle2D
    mode_indices: tuple[ModeIndex2D, ...]
    bloch_wave_vector: VectorQuantity

    def __post_init__(self) -> None:
        """Check every spectrum-request argument."""
        self._check_arg_model()
        self._check_arg_mode_indices()
        self._check_arg_bloch_wave_vector()

    def _check_arg_model(self) -> None:
        if not isinstance(self.model, PeriodicFreeParticle2D):
            raise TypeError("model must be PeriodicFreeParticle2D")

    def _check_arg_mode_indices(self) -> None:
        if type(self.mode_indices) is not tuple:
            raise TypeError("mode_indices must be a tuple")
        if len(self.mode_indices) == 0:
            raise ValueError("mode_indices must not be empty")
        for mode in self.mode_indices:
            if type(mode) is not tuple or len(mode) != 2:
                raise TypeError("each mode index must be a two-integer tuple")
            if any(type(component) is not int for component in mode):
                raise TypeError("mode-index components must be built-in integers")
        if len(set(self.mode_indices)) != len(self.mode_indices):
            raise ValueError("mode_indices must be unique")

    def _check_arg_bloch_wave_vector(self) -> None:
        if not isinstance(self.bloch_wave_vector, VectorQuantity):
            raise TypeError("bloch_wave_vector must be VectorQuantity")
        expected_unit = self.model.cell.reciprocal_basis.unit
        if self.bloch_wave_vector.unit != expected_unit:
            raise ValueError(
                "bloch_wave_vector must use the inverse selected length unit"
            )
        if self.bloch_wave_vector.magnitude.shape != (2,):
            raise ValueError("bloch_wave_vector must have shape (2,)")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicFreeParticle2DSpectrum(ResultsObject):
    """Return wave vectors and energies in deterministic request order."""

    request: PeriodicFreeParticle2DSpectrumRequest
    wave_vectors: MatrixQuantity
    energies: VectorQuantity

    def __post_init__(self) -> None:
        """Check every spectrum-result argument."""
        self._check_arg_request()
        self._check_arg_wave_vectors()
        self._check_arg_energies()

    def _check_arg_request(self) -> None:
        if not isinstance(self.request, PeriodicFreeParticle2DSpectrumRequest):
            raise TypeError("request must be PeriodicFreeParticle2DSpectrumRequest")

    def _check_arg_wave_vectors(self) -> None:
        if not isinstance(self.wave_vectors, MatrixQuantity):
            raise TypeError("wave_vectors must be MatrixQuantity")
        if self.wave_vectors.unit != self.request.model.cell.reciprocal_basis.unit:
            raise ValueError("wave_vectors must use the model reciprocal-length unit")
        if self.wave_vectors.magnitude.shape != (len(self.request.mode_indices), 2):
            raise ValueError("wave_vectors shape must match the requested modes")

    def _check_arg_energies(self) -> None:
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if self.energies.unit != self.request.model.cell.unit_system.energy_unit:
            raise ValueError("energies must use the selected energy unit")
        if self.energies.magnitude.shape != (len(self.request.mode_indices),):
            raise ValueError("energies shape must match the requested modes")
        if np.any(self.energies.magnitude < 0.0):
            raise ValueError("energies must be nonnegative")


class PeriodicFreeParticle2DSpectrumEvaluator(
    DataObjectActionizer[
        PeriodicFreeParticle2DSpectrumRequest,
        PeriodicFreeParticle2DSpectrum,
    ]
):
    """Evaluate exact Bloch-mode kinetic energies in the selected unit system."""

    def action(
        self,
        *,
        request: PeriodicFreeParticle2DSpectrumRequest,
    ) -> PeriodicFreeParticle2DSpectrum:
        """Return wave vectors and effective-mass kinetic energies."""
        if not isinstance(request, PeriodicFreeParticle2DSpectrumRequest):
            raise TypeError("request must be PeriodicFreeParticle2DSpectrumRequest")
        model = request.model
        reciprocal_basis = model.cell.reciprocal_basis
        mode_matrix = np.asarray(request.mode_indices, dtype=np.float64)
        wave_vectors = (
            mode_matrix @ reciprocal_basis.magnitude.T
            + request.bloch_wave_vector.magnitude[np.newaxis, :]
        )
        inverse_mass = np.linalg.inv(model.effective_mass.tensor.magnitude)
        quadratic_forms = np.einsum(
            "ni,ij,nj->n",
            wave_vectors,
            inverse_mass,
            wave_vectors,
        )
        hbar = model.cell.unit_system.hbar
        represented_energy_unit = PhysicalUnit(
            expression=(
                f"({hbar.unit.expression}) ** 2 / "
                f"(({model.effective_mass.tensor.unit.expression}) * "
                f"({model.cell.primitive_basis.unit.expression}) ** 2)"
            )
        )
        conversion_factor = MODEL_SYSTEM_UNIT_CONVERTER.conversion_factor(
            represented_energy_unit,
            model.cell.unit_system.energy_unit,
        )
        energies = 0.5 * hbar.magnitude**2 * quadratic_forms * conversion_factor
        return PeriodicFreeParticle2DSpectrum(
            request=request,
            wave_vectors=MatrixQuantity(
                magnitude=np.asarray(wave_vectors, dtype=np.float64),
                unit=reciprocal_basis.unit,
            ),
            energies=VectorQuantity(
                magnitude=np.asarray(energies, dtype=np.float64),
                unit=model.cell.unit_system.energy_unit,
            ),
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicFreeParticle2D(
    ActionizedDataObject[
        PeriodicFreeParticle2DSpectrumRequest,
        PeriodicFreeParticle2DSpectrum,
    ]
):
    """Define a free particle on a periodic two-dimensional primitive cell."""

    cell: PeriodicPrimitiveCell2D
    effective_mass: CartesianEffectiveMassTensor2D
    actionizer: DataObjectActionizer[
        PeriodicFreeParticle2DSpectrumRequest,
        PeriodicFreeParticle2DSpectrum,
    ] = field(
        default_factory=PeriodicFreeParticle2DSpectrumEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check every free-particle model argument."""
        self._check_arg_cell()
        self._check_arg_effective_mass()

    def _check_arg_cell(self) -> None:
        if not isinstance(self.cell, PeriodicPrimitiveCell2D):
            raise TypeError("cell must be PeriodicPrimitiveCell2D")

    def _check_arg_effective_mass(self) -> None:
        if not isinstance(self.effective_mass, CartesianEffectiveMassTensor2D):
            raise TypeError("effective_mass must be CartesianEffectiveMassTensor2D")
        if self.effective_mass.unit_system is not self.cell.unit_system:
            raise ValueError("cell and effective_mass must use the same unit system")

    def evaluate_modes(
        self,
        *,
        mode_indices: tuple[ModeIndex2D, ...],
        bloch_wave_vector: VectorQuantity,
    ) -> PeriodicFreeParticle2DSpectrum:
        """Evaluate requested reciprocal modes at one Bloch wave vector."""
        return self._respond(
            request=PeriodicFreeParticle2DSpectrumRequest(
                model=self,
                mode_indices=mode_indices,
                bloch_wave_vector=bloch_wave_vector,
            )
        )
