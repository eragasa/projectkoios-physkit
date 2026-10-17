"""Morse sampled evaluator."""

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.solidstate.interatomic_potentials.pair.base import (
    RadialPairPotentialEvaluation,
    RadialPairPotentialEvaluationRequest,
)
from projectkoios.physkit.solidstate.interatomic_potentials.pair.morse.potential import (  # noqa: E501
    MorsePotential,
)
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True)
class MorsePotentialEvaluator:
    """Evaluate Morse energies and $F_r=-dV/dr$."""

    def action(
        self,
        *,
        request: RadialPairPotentialEvaluationRequest,
    ) -> RadialPairPotentialEvaluation:
        """Evaluate in parameter-native units and convert the outputs."""
        if not isinstance(request, RadialPairPotentialEvaluationRequest):
            raise TypeError("request must be RadialPairPotentialEvaluationRequest")
        model = request.model
        if not isinstance(model, MorsePotential):
            raise TypeError("request model must be MorsePotential")
        parameters = model.parameters
        distance_unit = parameters.equilibrium_distance.unit
        distances = ModelSystemUnitConverter.convert_vector(
            request.distances,
            distance_unit,
        )
        inverse_distance_unit = PhysicalUnit(
            expression=f"1 / ({distance_unit.expression})"
        )
        inverse_range = ModelSystemUnitConverter.convert_scalar(
            parameters.inverse_range,
            inverse_distance_unit,
        )
        exponential = np.exp(
            -inverse_range.magnitude
            * (distances.magnitude - parameters.equilibrium_distance.magnitude)
        )
        energies = parameters.dissociation_energy.magnitude * (
            (1.0 - exponential) ** 2 - 1.0
        )
        forces = (
            -2.0
            * inverse_range.magnitude
            * parameters.dissociation_energy.magnitude
            * (1.0 - exponential)
            * exponential
        )
        native_force_unit = PhysicalUnit(
            expression=(
                f"({parameters.dissociation_energy.unit.expression}) / "
                f"({distance_unit.expression})"
            )
        )
        return RadialPairPotentialEvaluation(
            request=request,
            potential_energies=ModelSystemUnitConverter.convert_vector(
                VectorQuantity(
                    magnitude=np.asarray(energies, dtype=np.float64),
                    unit=parameters.dissociation_energy.unit,
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
