"""Analytical stationary eigenfunctions for the rectangular PIAB2D model."""

from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    MatrixQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    Unitless,
    VectorQuantity,
)

from ...base import Piab2D


@dataclass(frozen=True, slots=True, kw_only=True)
class Piab2DAnalyticalEigenfunctionEvaluationRequest(DataObject):
    """Request sampling of one analytical PIAB2D eigenvector."""

    eigenfunction: Piab2DAnalyticalEigenfunction
    coordinates_x: VectorQuantity
    coordinates_y: VectorQuantity

    def __post_init__(self) -> None:
        """Check each eigenfunction-evaluation request argument."""
        self._check_arg_eigenfunction()
        self._check_arg_coordinates_x()
        self._check_arg_coordinates_y()

    def _check_arg_eigenfunction(self) -> None:
        """Require the analytical eigenfunction that will be sampled."""
        if not isinstance(self.eigenfunction, Piab2DAnalyticalEigenfunction):
            raise TypeError("eigenfunction must be Piab2DAnalyticalEigenfunction")

    def _check_arg_coordinates_x(self) -> None:
        """Require vector-valued sampling coordinates along the x axis."""
        if not isinstance(self.coordinates_x, VectorQuantity):
            raise TypeError("coordinates_x must be VectorQuantity")

    def _check_arg_coordinates_y(self) -> None:
        """Require vector-valued sampling coordinates along the y axis."""
        if not isinstance(self.coordinates_y, VectorQuantity):
            raise TypeError("coordinates_y must be VectorQuantity")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class Piab2DAnalyticalEigenfunctionEvaluation(ResultsObject):
    """Correlate one analytical eigenvector with its sampled representation."""

    request: Piab2DAnalyticalEigenfunctionEvaluationRequest
    coordinates_x: VectorQuantity
    coordinates_y: VectorQuantity
    amplitudes: MatrixQuantity

    def __post_init__(self) -> None:
        """Check each sampled-eigenfunction response argument."""
        self._check_arg_request()
        self._check_arg_coordinates_x()
        self._check_arg_coordinates_y()
        self._check_arg_amplitudes()

    def _check_arg_request(self) -> None:
        """Require the complete analytical-eigenfunction evaluation request."""
        if not isinstance(self.request, Piab2DAnalyticalEigenfunctionEvaluationRequest):
            raise TypeError(
                "request must be Piab2DAnalyticalEigenfunctionEvaluationRequest"
            )

    def _check_arg_coordinates_x(self) -> None:
        """Require nonempty model-unit x coordinates inside the closed box."""
        if not isinstance(self.coordinates_x, VectorQuantity):
            raise TypeError("coordinates_x must be VectorQuantity")
        model = self.request.eigenfunction.model
        if self.coordinates_x.unit != model.unit_system.length_unit:
            raise ValueError("coordinates_x must use the model length unit")
        if self.coordinates_x.magnitude.size == 0:
            raise ValueError("coordinates_x must be nonempty")
        if np.any(self.coordinates_x.magnitude < 0.0) or np.any(
            self.coordinates_x.magnitude > model.length_x
        ):
            raise ValueError("coordinates_x must lie in the closed box interval")

    def _check_arg_coordinates_y(self) -> None:
        """Require nonempty model-unit y coordinates inside the closed box."""
        if not isinstance(self.coordinates_y, VectorQuantity):
            raise TypeError("coordinates_y must be VectorQuantity")
        model = self.request.eigenfunction.model
        if self.coordinates_y.unit != model.unit_system.length_unit:
            raise ValueError("coordinates_y must use the model length unit")
        if self.coordinates_y.magnitude.size == 0:
            raise ValueError("coordinates_y must be nonempty")
        if np.any(self.coordinates_y.magnitude < 0.0) or np.any(
            self.coordinates_y.magnitude > model.length_y
        ):
            raise ValueError("coordinates_y must lie in the closed box interval")

    def _check_arg_amplitudes(self) -> None:
        """Require normalized amplitudes on the complete rectangular grid."""
        if not isinstance(self.amplitudes, MatrixQuantity):
            raise TypeError("amplitudes must be MatrixQuantity")
        expected_shape = (
            self.coordinates_x.magnitude.size,
            self.coordinates_y.magnitude.size,
        )
        if self.amplitudes.magnitude.shape != expected_shape:
            raise ValueError("amplitudes must match the rectangular grid")
        length_unit = self.request.eigenfunction.model.unit_system.length_unit
        if isinstance(length_unit, Unitless):
            expected_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(length_unit, PhysicalUnit):
                raise TypeError("physical model requires a physical length unit")
            expected_unit = PhysicalUnit(expression=f"({length_unit.expression}) ** -1")
        if self.amplitudes.unit != expected_unit:
            raise ValueError("amplitudes must use the normalized eigenfunction unit")

    @property
    def eigenfunction(self) -> Piab2DAnalyticalEigenfunction:
        """Return the stationary eigenvector represented by this evaluation."""
        return self.request.eigenfunction

    @property
    def model(self) -> Piab2D:
        """Return the PIAB2D model that owns the eigenvector."""
        return self.eigenfunction.model

    @property
    def quantum_number_x(self) -> int:
        """Return the eigenvector's x quantum number."""
        return self.eigenfunction.quantum_number_x

    @property
    def quantum_number_y(self) -> int:
        """Return the eigenvector's y quantum number."""
        return self.eigenfunction.quantum_number_y

    @property
    def probability_density(self) -> MatrixQuantity:
        r"""Return the sampled probability density $|\phi(x,y)|^2$."""
        if isinstance(self.amplitudes.unit, Unitless):
            density_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(self.amplitudes.unit, PhysicalUnit):
                raise TypeError("amplitudes must use a physical or unitless unit")
            density_unit = PhysicalUnit(
                expression=f"({self.amplitudes.unit.expression}) ** 2"
            )
        return MatrixQuantity(
            magnitude=self.amplitudes.magnitude**2,
            unit=density_unit,
        )


@dataclass(frozen=True, slots=True)
class Piab2DAnalyticalEigenfunctionEvaluator:
    """Sample one normalized PIAB2D product eigenvector."""

    def action(
        self,
        *,
        request: Piab2DAnalyticalEigenfunctionEvaluationRequest,
    ) -> Piab2DAnalyticalEigenfunctionEvaluation:
        """Return the requested eigenvector on the rectangular coordinate grid."""
        if not isinstance(request, Piab2DAnalyticalEigenfunctionEvaluationRequest):
            raise TypeError(
                "request must be Piab2DAnalyticalEigenfunctionEvaluationRequest"
            )
        eigenfunction = request.eigenfunction
        model = eigenfunction.model
        selected_x = ModelSystemUnitConverter.convert_vector(
            request.coordinates_x,
            model.unit_system.length_unit,
        )
        selected_y = ModelSystemUnitConverter.convert_vector(
            request.coordinates_y,
            model.unit_system.length_unit,
        )
        if selected_x.magnitude.size == 0:
            raise ValueError("coordinates_x must be nonempty")
        if selected_y.magnitude.size == 0:
            raise ValueError("coordinates_y must be nonempty")
        if np.any(selected_x.magnitude < 0.0) or np.any(
            selected_x.magnitude > model.length_x
        ):
            raise ValueError("coordinates_x must lie in the closed box interval")
        if np.any(selected_y.magnitude < 0.0) or np.any(
            selected_y.magnitude > model.length_y
        ):
            raise ValueError("coordinates_y must lie in the closed box interval")

        x_factor = np.sin(
            eigenfunction.quantum_number_x
            * np.pi
            * selected_x.magnitude
            / model.length_x
        )
        y_factor = np.sin(
            eigenfunction.quantum_number_y
            * np.pi
            * selected_y.magnitude
            / model.length_y
        )
        values = (
            2.0
            / np.sqrt(model.length_x * model.length_y)
            * np.outer(x_factor, y_factor)
        )
        length_unit = model.unit_system.length_unit
        if isinstance(length_unit, Unitless):
            amplitude_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(length_unit, PhysicalUnit):
                raise TypeError("physical model requires a physical length unit")
            amplitude_unit = PhysicalUnit(
                expression=f"({length_unit.expression}) ** -1"
            )
        return Piab2DAnalyticalEigenfunctionEvaluation(
            request=request,
            coordinates_x=selected_x,
            coordinates_y=selected_y,
            amplitudes=MatrixQuantity(
                magnitude=values,
                unit=amplitude_unit,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class Piab2DAnalyticalEigenfunction(
    ActionizedDataObject[
        Piab2DAnalyticalEigenfunctionEvaluationRequest,
        Piab2DAnalyticalEigenfunctionEvaluation,
    ]
):
    """Represent one stationary PIAB2D Hamiltonian eigenvector."""

    model: Piab2D
    quantum_number_x: int
    quantum_number_y: int
    actionizer: Piab2DAnalyticalEigenfunctionEvaluator = field(
        default_factory=Piab2DAnalyticalEigenfunctionEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check each analytical-eigenfunction argument."""
        self._check_arg_model()
        self._check_arg_quantum_number_x()
        self._check_arg_quantum_number_y()

    def _check_arg_model(self) -> None:
        """Require the rectangular PIAB2D Hamiltonian model."""
        if not isinstance(self.model, Piab2D):
            raise TypeError("model must be Piab2D")

    def _check_arg_quantum_number_x(self) -> None:
        """Require a positive built-in integer x quantum number."""
        if type(self.quantum_number_x) is not int:
            raise TypeError("quantum_number_x must be a built-in integer")
        if self.quantum_number_x <= 0:
            raise ValueError("quantum_number_x must be positive")

    def _check_arg_quantum_number_y(self) -> None:
        """Require a positive built-in integer y quantum number."""
        if type(self.quantum_number_y) is not int:
            raise TypeError("quantum_number_y must be a built-in integer")
        if self.quantum_number_y <= 0:
            raise ValueError("quantum_number_y must be positive")

    def evaluate(
        self,
        *,
        coordinates_x: VectorQuantity,
        coordinates_y: VectorQuantity,
    ) -> Piab2DAnalyticalEigenfunctionEvaluation:
        """Sample this eigenvector on a rectangular coordinate grid."""
        return self._respond(
            request=Piab2DAnalyticalEigenfunctionEvaluationRequest(
                eigenfunction=self,
                coordinates_x=coordinates_x,
                coordinates_y=coordinates_y,
            )
        )
