"""Lennard--Jones 12-6 sampled evaluator."""

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.solidstate.interatomic_potentials.pair.base import (
    RadialPairPotentialEvaluation,
    RadialPairPotentialEvaluationRequest,
)
from projectkoios.physkit.solidstate.interatomic_potentials.pair.lennard_jones_12_6.potential import (  # noqa: E501
    LennardJonesPotential,
)
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True)
class LennardJonesPotentialEvaluator:
    """Evaluate Lennard--Jones energies and $F_r=-dV/dr$."""

    def action(
        self,
        *,
        request: RadialPairPotentialEvaluationRequest,
    ) -> RadialPairPotentialEvaluation:
        """Evaluate in parameter-native units and convert the outputs."""
        if not isinstance(request, RadialPairPotentialEvaluationRequest):
            raise TypeError("request must be RadialPairPotentialEvaluationRequest")
        model = request.model
        if not isinstance(model, LennardJonesPotential):
            raise TypeError("request model must be LennardJonesPotential")
        parameters = model.parameters
        distance_unit = parameters.zero_crossing_distance.unit
        distances = ModelSystemUnitConverter.convert_vector(
            request.distances,
            distance_unit,
        )
        ratio = parameters.zero_crossing_distance.magnitude / distances.magnitude
        ratio_sixth = ratio**6
        energies = (
            4.0 * parameters.well_depth.magnitude * (ratio_sixth**2 - ratio_sixth)
        )
        forces = (
            24.0
            * parameters.well_depth.magnitude
            / distances.magnitude
            * (2.0 * ratio_sixth**2 - ratio_sixth)
        )
        native_force_unit = PhysicalUnit(
            expression=(
                f"({parameters.well_depth.unit.expression}) / "
                f"({distance_unit.expression})"
            )
        )
        return RadialPairPotentialEvaluation(
            request=request,
            potential_energies=ModelSystemUnitConverter.convert_vector(
                VectorQuantity(
                    magnitude=np.asarray(energies, dtype=np.float64),
                    unit=parameters.well_depth.unit,
                ),
                request.energy_unit,
            ),
            radial_forces=ModelSystemUnitConverter.convert_vector(
                VectorQuantity(
                    magnitude=np.asarray(forces, dtype=np.float64),
                    unit=native_force_unit,
                ),
                request.force_unit,
            ),
        )
