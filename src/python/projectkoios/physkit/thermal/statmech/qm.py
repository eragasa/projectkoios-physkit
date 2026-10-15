r"""Finite-level canonical and Fermi--Dirac statistical analysis.

The actions in this module operate on explicit energy quantities. Thermal
energy, chemical potential, and represented level energies must use the same
unit. Weighted inversion additionally requires one explicit positive state
weight per level and a target in that weight unit. A Fermi--Dirac occupation
refers to one represented single-particle state; spin or other degeneracy is not
applied implicitly.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.qm.eigenfunctions import SampledEigenfunctions1D
from projectkoios.physkit.units import (
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, kw_only=True)
class FermiDiracOccupationEvaluationRequest(DataObject):
    """Request occupations for explicit energies under one thermal state."""

    state: FermiDiracState
    energies: VectorQuantity

    def __post_init__(self) -> None:
        """Check each occupation-evaluation request argument."""
        self._check_arg_state()
        self._check_arg_energies()

    def _check_arg_state(self) -> None:
        """Require the Fermi--Dirac state defining temperature and potential."""
        if not isinstance(self.state, FermiDiracState):
            raise TypeError("state must be FermiDiracState")

    def _check_arg_energies(self) -> None:
        """Require represented energies in the state's energy unit."""
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if self.energies.unit != self.state.thermal_energy.unit:
            raise ValueError("energies must use the Fermi-Dirac energy unit")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class FermiDiracOccupationEvaluation(ResultsObject):
    """Correlate an occupation request with represented occupations."""

    request: FermiDiracOccupationEvaluationRequest
    occupations: VectorQuantity

    def __post_init__(self) -> None:
        """Check each occupation-evaluation response argument."""
        self._check_arg_request()
        self._check_arg_occupations()

    def _check_arg_request(self) -> None:
        """Require the complete occupation-evaluation request."""
        if not isinstance(self.request, FermiDiracOccupationEvaluationRequest):
            raise TypeError("request must be FermiDiracOccupationEvaluationRequest")

    def _check_arg_occupations(self) -> None:
        """Require one unitless bounded occupation per requested energy."""
        if not isinstance(self.occupations, VectorQuantity):
            raise TypeError("occupations must be VectorQuantity")
        if not isinstance(self.occupations.unit, Unitless):
            raise ValueError("occupations must be unitless")
        if self.occupations.magnitude.shape != self.request.energies.magnitude.shape:
            raise ValueError("occupations must match energies")
        if np.any(self.occupations.magnitude < 0.0) or np.any(
            self.occupations.magnitude > 1.0
        ):
            raise ValueError("occupations must lie in the closed interval [0, 1]")


@dataclass(frozen=True, slots=True)
class FermiDiracOccupationAction:
    r"""Evaluate $[\exp((E-\mu)/(k_B T))+1]^{-1}$ stably."""

    def action(
        self,
        *,
        request: FermiDiracOccupationEvaluationRequest,
    ) -> FermiDiracOccupationEvaluation:
        """Return one occupation for every represented energy level."""
        if not isinstance(request, FermiDiracOccupationEvaluationRequest):
            raise TypeError("request must be FermiDiracOccupationEvaluationRequest")
        arguments = (
            request.energies.magnitude - request.state.chemical_potential.magnitude
        ) / request.state.thermal_energy.magnitude
        occupations = np.empty(arguments.shape, dtype=np.float64)
        nonnegative = arguments >= 0.0
        decaying = np.exp(-arguments[nonnegative])
        occupations[nonnegative] = decaying / (1.0 + decaying)
        increasing = np.exp(arguments[~nonnegative])
        occupations[~nonnegative] = 1.0 / (1.0 + increasing)
        return FermiDiracOccupationEvaluation(
            request=request,
            occupations=VectorQuantity(
                magnitude=occupations,
                unit=Unitless(),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class FermiDiracState(
    ActionizedDataObject[
        FermiDiracOccupationEvaluationRequest,
        FermiDiracOccupationEvaluation,
    ]
):
    """Represent thermal energy and chemical potential in one energy unit."""

    thermal_energy: ScalarQuantity
    chemical_potential: ScalarQuantity
    actionizer: FermiDiracOccupationAction = field(
        default_factory=FermiDiracOccupationAction,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check each Fermi--Dirac state argument."""
        self._check_arg_thermal_energy()
        self._check_arg_chemical_potential()

    def _check_arg_thermal_energy(self) -> None:
        """Require a positive scalar thermal energy."""
        if not isinstance(self.thermal_energy, ScalarQuantity):
            raise TypeError("thermal_energy must be ScalarQuantity")
        if self.thermal_energy.magnitude <= 0.0:
            raise ValueError("thermal_energy must be positive")

    def _check_arg_chemical_potential(self) -> None:
        """Require a scalar chemical potential in the thermal-energy unit."""
        if not isinstance(self.chemical_potential, ScalarQuantity):
            raise TypeError("chemical_potential must be ScalarQuantity")
        if self.chemical_potential.unit != self.thermal_energy.unit:
            raise ValueError(
                "chemical_potential and thermal_energy must use the same unit"
            )

    def evaluate_occupations(
        self,
        *,
        energies: VectorQuantity,
    ) -> FermiDiracOccupationEvaluation:
        """Evaluate occupations for explicit represented energies."""
        return self._respond(
            request=FermiDiracOccupationEvaluationRequest(
                state=self,
                energies=energies,
            )
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class FermiDiracChemicalPotentialSolveRequest(DataObject):
    """Request a finite-level chemical potential satisfying a target."""

    solver: FermiDiracChemicalPotentialSolver
    energies: VectorQuantity
    thermal_energy: ScalarQuantity
    target_occupation: ScalarQuantity

    def __post_init__(self) -> None:
        """Check each chemical-potential solve-request argument."""
        self._check_arg_solver()
        self._check_arg_energies()
        self._check_arg_thermal_energy()
        self._check_arg_target_occupation()

    def _check_arg_solver(self) -> None:
        """Require the configured bisection solver façade."""
        if not isinstance(self.solver, FermiDiracChemicalPotentialSolver):
            raise TypeError("solver must be FermiDiracChemicalPotentialSolver")

    def _check_arg_energies(self) -> None:
        """Require a nonempty finite-level energy inventory."""
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if self.energies.magnitude.size == 0:
            raise ValueError("energies must be nonempty")

    def _check_arg_thermal_energy(self) -> None:
        """Require positive thermal energy in the represented energy unit."""
        if not isinstance(self.thermal_energy, ScalarQuantity):
            raise TypeError("thermal_energy must be ScalarQuantity")
        if self.thermal_energy.magnitude <= 0.0:
            raise ValueError("thermal_energy must be positive")
        if self.thermal_energy.unit != self.energies.unit:
            raise ValueError("thermal_energy and energies must use the same unit")

    def _check_arg_target_occupation(self) -> None:
        """Require a unitless target strictly between zero and level count."""
        if not isinstance(self.target_occupation, ScalarQuantity):
            raise TypeError("target_occupation must be ScalarQuantity")
        if not isinstance(self.target_occupation.unit, Unitless):
            raise ValueError("target_occupation must be unitless")
        if (
            not 0.0
            < self.target_occupation.magnitude
            < float(self.energies.magnitude.size)
        ):
            raise ValueError(
                "target_occupation must lie strictly between zero and level count"
            )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class FermiDiracChemicalPotentialResult(ResultsObject):
    """Correlate a solve request with chemical potential and occupations."""

    request: FermiDiracChemicalPotentialSolveRequest
    chemical_potential: ScalarQuantity
    occupations: VectorQuantity
    occupation_residual: ScalarQuantity

    def __post_init__(self) -> None:
        """Check each chemical-potential response argument."""
        self._check_arg_request()
        self._check_arg_chemical_potential()
        self._check_arg_occupations()
        self._check_arg_occupation_residual()

    def _check_arg_request(self) -> None:
        """Require the complete chemical-potential solve request."""
        if not isinstance(self.request, FermiDiracChemicalPotentialSolveRequest):
            raise TypeError("request must be FermiDiracChemicalPotentialSolveRequest")

    def _check_arg_chemical_potential(self) -> None:
        """Require the solved potential in the represented energy unit."""
        if not isinstance(self.chemical_potential, ScalarQuantity):
            raise TypeError("chemical_potential must be ScalarQuantity")
        if self.chemical_potential.unit != self.request.energies.unit:
            raise ValueError("chemical_potential must use the represented energy unit")

    def _check_arg_occupations(self) -> None:
        """Require one unitless bounded occupation per represented level."""
        if not isinstance(self.occupations, VectorQuantity):
            raise TypeError("occupations must be VectorQuantity")
        if not isinstance(self.occupations.unit, Unitless):
            raise ValueError("occupations must be unitless")
        if self.occupations.magnitude.shape != self.request.energies.magnitude.shape:
            raise ValueError("occupations must match energies")
        if np.any(self.occupations.magnitude < 0.0) or np.any(
            self.occupations.magnitude > 1.0
        ):
            raise ValueError("occupations must lie in the closed interval [0, 1]")

    def _check_arg_occupation_residual(self) -> None:
        """Require the unitless residual implied by occupations and target."""
        if not isinstance(self.occupation_residual, ScalarQuantity):
            raise TypeError("occupation_residual must be ScalarQuantity")
        if not isinstance(self.occupation_residual.unit, Unitless):
            raise ValueError("occupation_residual must be unitless")
        expected_residual = float(
            np.sum(self.occupations.magnitude)
            - self.request.target_occupation.magnitude
        )
        tolerance = 16.0 * np.finfo(np.float64).eps * max(1.0, abs(expected_residual))
        if not np.isclose(
            self.occupation_residual.magnitude,
            expected_residual,
            rtol=0.0,
            atol=tolerance,
        ):
            raise ValueError(
                "occupation_residual must match the represented occupations"
            )

    @property
    def target_occupation(self) -> ScalarQuantity:
        """Return the requested unitless total occupation."""
        return self.request.target_occupation


@dataclass(frozen=True, slots=True, kw_only=True)
class WeightedFermiDiracChemicalPotentialSolveRequest(DataObject):
    """Request a chemical potential for explicitly weighted states.

    Parameters
    ----------
    solver
        Existing deterministic Fermi--Dirac chemical-potential solver whose
        tolerances and iteration limit control the weighted solve.
    energies
        One represented single-particle energy per state.
    state_weights
        Positive physical or unitless integration weight for every energy.
        Spin, valley, quadrature, and density normalization must already be
        represented here.
    thermal_energy
        Positive :math:`k_B T` in the exact energy unit.
    target_state_weight
        Requested occupied state weight in the state-weight unit.
    """

    solver: FermiDiracChemicalPotentialSolver
    energies: VectorQuantity
    state_weights: VectorQuantity
    thermal_energy: ScalarQuantity
    target_state_weight: ScalarQuantity

    def __post_init__(self) -> None:
        """Check the weighted state inventory and target."""
        if not isinstance(self.solver, FermiDiracChemicalPotentialSolver):
            raise TypeError("solver must be FermiDiracChemicalPotentialSolver")
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if self.energies.magnitude.size == 0:
            raise ValueError("energies must be nonempty")
        if not isinstance(self.state_weights, VectorQuantity):
            raise TypeError("state_weights must be VectorQuantity")
        if self.state_weights.magnitude.shape != self.energies.magnitude.shape:
            raise ValueError("state_weights must match energies")
        if not isinstance(self.state_weights.unit, (PhysicalUnit, Unitless)):
            raise TypeError("state_weights must use PhysicalUnit or Unitless")
        if np.any(self.state_weights.magnitude <= 0.0):
            raise ValueError("state_weights must be positive")
        if not isinstance(self.thermal_energy, ScalarQuantity):
            raise TypeError("thermal_energy must be ScalarQuantity")
        if self.thermal_energy.unit != self.energies.unit:
            raise ValueError("thermal_energy and energies must use the same unit")
        if self.thermal_energy.magnitude <= 0.0:
            raise ValueError("thermal_energy must be positive")
        if not isinstance(self.target_state_weight, ScalarQuantity):
            raise TypeError("target_state_weight must be ScalarQuantity")
        if self.target_state_weight.unit != self.state_weights.unit:
            raise ValueError("target_state_weight must use the state-weight unit")
        total_state_weight = float(np.sum(self.state_weights.magnitude))
        if not 0.0 < self.target_state_weight.magnitude < total_state_weight:
            raise ValueError(
                "target_state_weight must lie strictly between zero and total "
                "state weight"
            )

    @property
    def total_state_weight(self) -> ScalarQuantity:
        """Return the available weight represented by all states."""
        return ScalarQuantity(
            magnitude=float(np.sum(self.state_weights.magnitude)),
            unit=self.state_weights.unit,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class WeightedFermiDiracChemicalPotentialResult(ResultsObject):
    """Correlate a weighted solve with occupations and its residual.

    Parameters
    ----------
    request
        Complete weighted inversion request.
    chemical_potential
        Solved chemical potential in the represented energy unit.
    occupations
        Unitless Fermi--Dirac occupation for every represented state.
    occupied_state_weight
        Direct weighted sum of the returned occupations.
    state_weight_residual
        Numerically stable occupied-target residual in the state-weight unit.
    """

    request: WeightedFermiDiracChemicalPotentialSolveRequest
    chemical_potential: ScalarQuantity
    occupations: VectorQuantity
    occupied_state_weight: ScalarQuantity
    state_weight_residual: ScalarQuantity

    def __post_init__(self) -> None:
        """Check weighted chemical-potential result invariants."""
        if not isinstance(
            self.request,
            WeightedFermiDiracChemicalPotentialSolveRequest,
        ):
            raise TypeError(
                "request must be WeightedFermiDiracChemicalPotentialSolveRequest"
            )
        if not isinstance(self.chemical_potential, ScalarQuantity):
            raise TypeError("chemical_potential must be ScalarQuantity")
        if self.chemical_potential.unit != self.request.energies.unit:
            raise ValueError("chemical_potential must use the represented energy unit")
        if not isinstance(self.occupations, VectorQuantity):
            raise TypeError("occupations must be VectorQuantity")
        if not isinstance(self.occupations.unit, Unitless):
            raise ValueError("occupations must be unitless")
        if self.occupations.magnitude.shape != self.request.energies.magnitude.shape:
            raise ValueError("occupations must match energies")
        if np.any(self.occupations.magnitude < 0.0) or np.any(
            self.occupations.magnitude > 1.0
        ):
            raise ValueError("occupations must lie in the closed interval [0, 1]")
        expected_occupied = float(
            np.dot(
                self.request.state_weights.magnitude,
                self.occupations.magnitude,
            )
        )
        self._check_weight_scalar(
            value=self.occupied_state_weight,
            name="occupied_state_weight",
            expected_magnitude=expected_occupied,
        )
        total_weight = self.request.total_state_weight.magnitude
        target_weight = self.request.target_state_weight.magnitude
        # Near full filling, form the residual from holes rather than
        # subtracting two nearly equal extensive occupied-state densities.
        if target_weight <= 0.5 * total_weight:
            expected_residual = expected_occupied - target_weight
        else:
            expected_unoccupied = float(
                np.dot(
                    self.request.state_weights.magnitude,
                    1.0 - self.occupations.magnitude,
                )
            )
            expected_residual = total_weight - target_weight - expected_unoccupied
        self._check_weight_scalar(
            value=self.state_weight_residual,
            name="state_weight_residual",
            expected_magnitude=expected_residual,
        )

    def _check_weight_scalar(
        self,
        *,
        value: ScalarQuantity,
        name: str,
        expected_magnitude: float,
    ) -> None:
        """Check one state-weight result against the represented sum."""
        if not isinstance(value, ScalarQuantity):
            raise TypeError(f"{name} must be ScalarQuantity")
        if value.unit != self.request.state_weights.unit:
            raise ValueError(f"{name} must use the state-weight unit")
        tolerance = (
            32.0
            * np.finfo(np.float64).eps
            * max(
                1.0,
                abs(expected_magnitude),
            )
        )
        if not np.isclose(
            value.magnitude,
            expected_magnitude,
            rtol=0.0,
            atol=tolerance,
        ):
            raise ValueError(f"{name} must match the weighted occupations")

    @property
    def target_state_weight(self) -> ScalarQuantity:
        """Return the requested occupied state weight."""
        return self.request.target_state_weight


@dataclass(frozen=True, slots=True)
class FermiDiracChemicalPotentialBisection:
    """Solve the finite-level particle-number equation by bisection."""

    BRACKET_THERMAL_ENERGY_MULTIPLIER = 64.0

    def action(
        self,
        *,
        request: FermiDiracChemicalPotentialSolveRequest,
    ) -> FermiDiracChemicalPotentialResult:
        """Return the chemical potential matching the requested target."""
        if not isinstance(request, FermiDiracChemicalPotentialSolveRequest):
            raise TypeError("request must be FermiDiracChemicalPotentialSolveRequest")
        thermal_scale = request.thermal_energy.magnitude
        lower = float(
            np.min(request.energies.magnitude)
            - self.BRACKET_THERMAL_ENERGY_MULTIPLIER * thermal_scale
        )
        upper = float(
            np.max(request.energies.magnitude)
            + self.BRACKET_THERMAL_ENERGY_MULTIPLIER * thermal_scale
        )
        energy_unit = request.energies.unit
        midpoint = 0.5 * (lower + upper)
        occupations = VectorQuantity(
            magnitude=np.empty(request.energies.magnitude.shape, dtype=np.float64),
            unit=Unitless(),
        )
        residual = float("inf")

        for _ in range(request.solver.maximum_iterations):
            midpoint = 0.5 * (lower + upper)
            state = FermiDiracState(
                thermal_energy=request.thermal_energy,
                chemical_potential=ScalarQuantity(
                    magnitude=midpoint,
                    unit=energy_unit,
                ),
            )
            occupations = state.evaluate_occupations(
                energies=request.energies
            ).occupations
            residual = float(
                np.sum(occupations.magnitude) - request.target_occupation.magnitude
            )
            if abs(residual) <= request.solver.occupation_tolerance:
                break
            if residual < 0.0:
                lower = midpoint
            else:
                upper = midpoint
        else:
            raise RuntimeError(
                "chemical-potential bisection did not reach occupation_tolerance"
            )

        return FermiDiracChemicalPotentialResult(
            request=request,
            chemical_potential=ScalarQuantity(
                magnitude=midpoint,
                unit=energy_unit,
            ),
            occupations=occupations,
            occupation_residual=ScalarQuantity(
                magnitude=residual,
                unit=Unitless(),
            ),
        )


@dataclass(frozen=True, slots=True)
class WeightedFermiDiracChemicalPotentialBisection:
    """Solve a weighted Fermi--Dirac state inventory by bisection.

    Notes
    -----
    This action reuses :class:`FermiDiracState` for stable occupation
    evaluation. It expands the initial energy bracket when a target is very
    close to empty or full filling and evaluates residuals through the smaller
    particle or hole population to reduce cancellation.
    """

    BRACKET_THERMAL_ENERGY_MULTIPLIER = 64.0
    MAXIMUM_BRACKET_EXPANSIONS = 16

    def action(
        self,
        *,
        request: WeightedFermiDiracChemicalPotentialSolveRequest,
    ) -> WeightedFermiDiracChemicalPotentialResult:
        """Return the chemical potential matching a weighted target.

        Parameters
        ----------
        request
            Checked energy inventory, state weights, thermal energy, target,
            and deterministic solver controls.

        Returns
        -------
        WeightedFermiDiracChemicalPotentialResult
            Solved potential, occupations, occupied weight, and residual.
        """
        if not isinstance(
            request,
            WeightedFermiDiracChemicalPotentialSolveRequest,
        ):
            raise TypeError(
                "request must be WeightedFermiDiracChemicalPotentialSolveRequest"
            )
        energy_unit = request.energies.unit
        thermal_scale = request.thermal_energy.magnitude
        span = self.BRACKET_THERMAL_ENERGY_MULTIPLIER * thermal_scale
        minimum_energy = float(np.min(request.energies.magnitude))
        maximum_energy = float(np.max(request.energies.magnitude))
        total_weight = request.total_state_weight.magnitude
        target_weight = request.target_state_weight.magnitude

        def evaluate(mu: float) -> tuple[VectorQuantity, float, float]:
            state = FermiDiracState(
                thermal_energy=request.thermal_energy,
                chemical_potential=ScalarQuantity(
                    magnitude=mu,
                    unit=energy_unit,
                ),
            )
            occupations = state.evaluate_occupations(
                energies=request.energies
            ).occupations
            occupied = float(
                np.dot(
                    request.state_weights.magnitude,
                    occupations.magnitude,
                )
            )
            # Use particles below half filling and holes above half filling.
            # This avoids cancellation between total and occupied weights near
            # a nearly full inventory while preserving the residual sign.
            if target_weight <= 0.5 * total_weight:
                residual = occupied - target_weight
            else:
                unoccupied = float(
                    np.dot(
                        request.state_weights.magnitude,
                        1.0 - occupations.magnitude,
                    )
                )
                residual = total_weight - target_weight - unoccupied
            return occupations, occupied, residual

        # Sixty-four thermal energies normally saturate both bracket ends.
        # Expansion keeps extreme but attainable filling targets supported.
        lower = minimum_energy - span
        upper = maximum_energy + span
        _, _, lower_residual = evaluate(lower)
        _, _, upper_residual = evaluate(upper)
        for _ in range(self.MAXIMUM_BRACKET_EXPANSIONS):
            if lower_residual < 0.0 and upper_residual > 0.0:
                break
            span *= 2.0
            if lower_residual >= 0.0:
                lower = minimum_energy - span
                _, _, lower_residual = evaluate(lower)
            if upper_residual <= 0.0:
                upper = maximum_energy + span
                _, _, upper_residual = evaluate(upper)
        if lower_residual >= 0.0 or upper_residual <= 0.0:
            raise RuntimeError("weighted chemical-potential root was not bracketed")

        # Scale convergence by the smaller particle/hole population. A
        # binary64 floor prevents demanding resolution below roundoff.
        limiting_population = min(target_weight, total_weight - target_weight)
        residual_tolerance = (
            max(
                request.solver.occupation_tolerance,
                64.0 * np.finfo(np.float64).eps,
            )
            * limiting_population
        )
        midpoint = 0.5 * (lower + upper)
        occupations, occupied, residual = evaluate(midpoint)
        for _ in range(request.solver.maximum_iterations):
            midpoint = 0.5 * (lower + upper)
            occupations, occupied, residual = evaluate(midpoint)
            if abs(residual) <= residual_tolerance:
                break
            if residual < 0.0:
                lower = midpoint
            else:
                upper = midpoint
        else:
            raise RuntimeError(
                "weighted chemical-potential bisection did not reach the relative "
                "state-weight tolerance"
            )

        return WeightedFermiDiracChemicalPotentialResult(
            request=request,
            chemical_potential=ScalarQuantity(
                magnitude=midpoint,
                unit=energy_unit,
            ),
            occupations=occupations,
            occupied_state_weight=ScalarQuantity(
                magnitude=occupied,
                unit=request.state_weights.unit,
            ),
            state_weight_residual=ScalarQuantity(
                magnitude=residual,
                unit=request.state_weights.unit,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class FermiDiracChemicalPotentialSolver(
    ActionizedDataObject[
        FermiDiracChemicalPotentialSolveRequest,
        FermiDiracChemicalPotentialResult,
    ]
):
    """Configure deterministic finite-level chemical-potential bisection."""

    occupation_tolerance: float = 1.0e-12
    maximum_iterations: int = 256
    actionizer: FermiDiracChemicalPotentialBisection = field(
        default_factory=FermiDiracChemicalPotentialBisection,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check each deterministic bisection-control argument."""
        self._check_arg_occupation_tolerance()
        self._check_arg_maximum_iterations()

    def _check_arg_occupation_tolerance(self) -> None:
        """Require a finite positive built-in floating-point tolerance."""
        if type(self.occupation_tolerance) is not float:
            raise TypeError("occupation_tolerance must be a built-in float")
        if not np.isfinite(self.occupation_tolerance):
            raise ValueError("occupation_tolerance must be finite")
        if self.occupation_tolerance <= 0.0:
            raise ValueError("occupation_tolerance must be positive")

    def _check_arg_maximum_iterations(self) -> None:
        """Require a positive built-in integer iteration limit."""
        if type(self.maximum_iterations) is not int:
            raise TypeError("maximum_iterations must be a built-in int")
        if self.maximum_iterations <= 0:
            raise ValueError("maximum_iterations must be positive")

    def solve(
        self,
        *,
        energies: VectorQuantity,
        thermal_energy: ScalarQuantity,
        target_occupation: ScalarQuantity,
    ) -> FermiDiracChemicalPotentialResult:
        """Solve for the chemical potential matching an occupation target."""
        return self._respond(
            request=FermiDiracChemicalPotentialSolveRequest(
                solver=self,
                energies=energies,
                thermal_energy=thermal_energy,
                target_occupation=target_occupation,
            )
        )

    def solve_weighted(
        self,
        *,
        energies: VectorQuantity,
        state_weights: VectorQuantity,
        thermal_energy: ScalarQuantity,
        target_state_weight: ScalarQuantity,
    ) -> WeightedFermiDiracChemicalPotentialResult:
        """Solve for a chemical potential using explicit state weights.

        Parameters
        ----------
        energies
            One energy per represented single-particle state.
        state_weights
            Positive integration or density weight for each energy.
        thermal_energy
            Positive :math:`k_B T` in the energy unit.
        target_state_weight
            Requested occupied weight in the state-weight unit.

        Returns
        -------
        WeightedFermiDiracChemicalPotentialResult
            Unit-aware weighted inversion result.

        Notes
        -----
        Degeneracy is not introduced implicitly. Callers include spin, valley,
        or other multiplicities in ``state_weights``.
        """
        return WeightedFermiDiracChemicalPotentialBisection().action(
            request=WeightedFermiDiracChemicalPotentialSolveRequest(
                solver=self,
                energies=energies,
                state_weights=state_weights,
                thermal_energy=thermal_energy,
                target_state_weight=target_state_weight,
            )
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class CanonicalLevelProbabilityEvaluationRequest(DataObject):
    """Request canonical probabilities for a finite level inventory."""

    model: CanonicalLevelProbabilityModel
    energies: VectorQuantity

    def __post_init__(self) -> None:
        """Check each canonical-probability request argument."""
        self._check_arg_model()
        self._check_arg_energies()

    def _check_arg_model(self) -> None:
        """Require the canonical finite-level model."""
        if not isinstance(self.model, CanonicalLevelProbabilityModel):
            raise TypeError("model must be CanonicalLevelProbabilityModel")

    def _check_arg_energies(self) -> None:
        """Require nonempty energies in the model's thermal-energy unit."""
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if self.energies.magnitude.size == 0:
            raise ValueError("energies must be nonempty")
        if self.energies.unit != self.model.thermal_energy.unit:
            raise ValueError("energies must use the thermal-energy unit")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CanonicalLevelProbabilityEvaluation(ResultsObject):
    """Correlate a canonical request with normalized probabilities."""

    request: CanonicalLevelProbabilityEvaluationRequest
    probabilities: VectorQuantity

    def __post_init__(self) -> None:
        """Check each canonical-probability response argument."""
        self._check_arg_request()
        self._check_arg_probabilities()

    def _check_arg_request(self) -> None:
        """Require the complete canonical-probability request."""
        if not isinstance(self.request, CanonicalLevelProbabilityEvaluationRequest):
            raise TypeError(
                "request must be CanonicalLevelProbabilityEvaluationRequest"
            )

    def _check_arg_probabilities(self) -> None:
        """Require normalized unitless probabilities matching the levels."""
        if not isinstance(self.probabilities, VectorQuantity):
            raise TypeError("probabilities must be VectorQuantity")
        if not isinstance(self.probabilities.unit, Unitless):
            raise ValueError("probabilities must be unitless")
        if self.probabilities.magnitude.shape != self.request.energies.magnitude.shape:
            raise ValueError("probabilities must match energies")
        tolerance = 32.0 * np.finfo(np.float64).eps
        if not np.isclose(
            np.sum(self.probabilities.magnitude),
            1.0,
            rtol=tolerance,
            atol=tolerance,
        ):
            raise ValueError("probabilities must sum to one")


@dataclass(frozen=True, slots=True)
class CanonicalLevelProbabilityAction:
    """Evaluate stable normalized finite-level canonical probabilities."""

    def action(
        self,
        *,
        request: CanonicalLevelProbabilityEvaluationRequest,
    ) -> CanonicalLevelProbabilityEvaluation:
        """Return normalized Boltzmann probabilities for the requested levels."""
        if not isinstance(request, CanonicalLevelProbabilityEvaluationRequest):
            raise TypeError(
                "request must be CanonicalLevelProbabilityEvaluationRequest"
            )
        shifted = request.energies.magnitude - float(np.min(request.energies.magnitude))
        weights = np.exp(-shifted / request.model.thermal_energy.magnitude)
        return CanonicalLevelProbabilityEvaluation(
            request=request,
            probabilities=VectorQuantity(
                magnitude=weights / np.sum(weights),
                unit=Unitless(),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class CanonicalLevelProbabilityModel(
    ActionizedDataObject[
        CanonicalLevelProbabilityEvaluationRequest,
        CanonicalLevelProbabilityEvaluation,
    ]
):
    """Represent the thermal energy for finite-level canonical probabilities."""

    thermal_energy: ScalarQuantity
    actionizer: CanonicalLevelProbabilityAction = field(
        default_factory=CanonicalLevelProbabilityAction,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the thermal-energy argument."""
        self._check_arg_thermal_energy()

    def _check_arg_thermal_energy(self) -> None:
        """Require a positive scalar thermal energy."""
        if not isinstance(self.thermal_energy, ScalarQuantity):
            raise TypeError("thermal_energy must be ScalarQuantity")
        if self.thermal_energy.magnitude <= 0.0:
            raise ValueError("thermal_energy must be positive")

    def evaluate(
        self,
        *,
        energies: VectorQuantity,
    ) -> CanonicalLevelProbabilityEvaluation:
        """Evaluate normalized probabilities for explicit energy levels."""
        return self._respond(
            request=CanonicalLevelProbabilityEvaluationRequest(
                model=self,
                energies=energies,
            )
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class OccupiedEigenfunctionDensityEvaluationRequest(DataObject):
    """Request an occupation-weighted density from sampled eigenfunctions."""

    eigenfunctions: SampledEigenfunctions1D
    occupations: VectorQuantity

    def __post_init__(self) -> None:
        """Check each occupied-density request argument."""
        self._check_arg_eigenfunctions()
        self._check_arg_occupations()

    def _check_arg_eigenfunctions(self) -> None:
        """Require an ordered sampled eigenfunction collection."""
        if not isinstance(self.eigenfunctions, SampledEigenfunctions1D):
            raise TypeError("eigenfunctions must be SampledEigenfunctions1D")

    def _check_arg_occupations(self) -> None:
        """Require one nonnegative unitless occupation per eigenfunction."""
        if not isinstance(self.occupations, VectorQuantity):
            raise TypeError("occupations must be VectorQuantity")
        if not isinstance(self.occupations.unit, Unitless):
            raise ValueError("occupations must be unitless")
        if self.occupations.magnitude.shape != (len(self.eigenfunctions),):
            raise ValueError("occupations must match the eigenfunctions")
        if np.any(self.occupations.magnitude < 0.0):
            raise ValueError("occupations must be nonnegative")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class OccupiedEigenfunctionDensityEvaluation(ResultsObject):
    """Correlate an occupied-density request with sampled density values."""

    request: OccupiedEigenfunctionDensityEvaluationRequest
    density: VectorQuantity

    def __post_init__(self) -> None:
        """Check each occupied-density response argument."""
        self._check_arg_request()
        self._check_arg_density()

    def _check_arg_request(self) -> None:
        """Require the complete occupied-density request."""
        if not isinstance(self.request, OccupiedEigenfunctionDensityEvaluationRequest):
            raise TypeError(
                "request must be OccupiedEigenfunctionDensityEvaluationRequest"
            )

    def _check_arg_density(self) -> None:
        """Require one nonnegative density value per sampled coordinate."""
        if not isinstance(self.density, VectorQuantity):
            raise TypeError("density must be VectorQuantity")
        coordinates = self.request.eigenfunctions.coordinates
        if self.density.magnitude.shape != coordinates.magnitude.shape:
            raise ValueError("density must match the eigenfunction coordinates")
        if np.any(self.density.magnitude < 0.0):
            raise ValueError("density must be nonnegative")


@dataclass(frozen=True, slots=True)
class OccupiedEigenfunctionDensityAction:
    r"""Evaluate $\rho(x_i)=\sum_n f_n|\phi_n(x_i)|^2$."""

    def action(
        self,
        *,
        request: OccupiedEigenfunctionDensityEvaluationRequest,
    ) -> OccupiedEigenfunctionDensityEvaluation:
        """Return occupation-weighted density on shared sample coordinates."""
        if not isinstance(request, OccupiedEigenfunctionDensityEvaluationRequest):
            raise TypeError(
                "request must be OccupiedEigenfunctionDensityEvaluationRequest"
            )
        amplitudes = request.eigenfunctions.amplitude_matrix
        density = np.abs(amplitudes.magnitude) ** 2 @ request.occupations.magnitude
        if isinstance(amplitudes.unit, Unitless):
            density_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(amplitudes.unit, PhysicalUnit):
                raise TypeError(
                    "eigenfunction amplitudes must use a physical or unitless unit"
                )
            density_unit = PhysicalUnit(
                expression=f"({amplitudes.unit.expression}) ** 2"
            )
        return OccupiedEigenfunctionDensityEvaluation(
            request=request,
            density=VectorQuantity(
                magnitude=density,
                unit=density_unit,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class OccupiedEigenfunctionDensityEvaluator(
    ActionizedDataObject[
        OccupiedEigenfunctionDensityEvaluationRequest,
        OccupiedEigenfunctionDensityEvaluation,
    ]
):
    """Evaluate occupation-weighted density for sampled eigenfunctions."""

    actionizer: OccupiedEigenfunctionDensityAction = field(
        default_factory=OccupiedEigenfunctionDensityAction,
        init=False,
        repr=False,
        compare=False,
    )

    def evaluate(
        self,
        *,
        eigenfunctions: SampledEigenfunctions1D,
        occupations: VectorQuantity,
    ) -> OccupiedEigenfunctionDensityEvaluation:
        """Evaluate density for explicit eigenfunctions and occupations."""
        return self._respond(
            request=OccupiedEigenfunctionDensityEvaluationRequest(
                eigenfunctions=eigenfunctions,
                occupations=occupations,
            )
        )
