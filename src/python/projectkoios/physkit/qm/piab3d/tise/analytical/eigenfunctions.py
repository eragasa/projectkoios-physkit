"""Analytical eigenfunction plane slices for the rectangular PIAB3D model."""

from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    MatrixQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

from ...base import Piab3D


class Piab3DPlaneNormal(Enum):
    """Identify the Cartesian axis held fixed by a PIAB3D plane slice."""

    X = "x"
    Y = "y"
    Z = "z"


@dataclass(frozen=True, slots=True, kw_only=True)
class Piab3DAnalyticalEigenfunctionPlaneSliceEvaluationRequest(DataObject):
    """Request one Cartesian plane slice through a PIAB3D eigenvector."""

    eigenfunction: Piab3DAnalyticalEigenfunction
    plane_normal: Piab3DPlaneNormal
    coordinates_u: VectorQuantity
    coordinates_v: VectorQuantity
    fixed_coordinate: ScalarQuantity

    def __post_init__(self) -> None:
        """Check each plane-slice evaluation request argument."""
        self._check_arg_eigenfunction()
        self._check_arg_plane_normal()
        self._check_arg_coordinates_u()
        self._check_arg_coordinates_v()
        self._check_arg_fixed_coordinate()

    def _check_arg_eigenfunction(self) -> None:
        """Require the stationary eigenvector whose slice will be sampled."""
        if not isinstance(self.eigenfunction, Piab3DAnalyticalEigenfunction):
            raise TypeError("eigenfunction must be Piab3DAnalyticalEigenfunction")

    def _check_arg_plane_normal(self) -> None:
        """Require one explicit Cartesian plane normal."""
        if not isinstance(self.plane_normal, Piab3DPlaneNormal):
            raise TypeError("plane_normal must be Piab3DPlaneNormal")

    def _check_arg_coordinates_u(self) -> None:
        """Require vector-valued coordinates along the first in-plane axis."""
        if not isinstance(self.coordinates_u, VectorQuantity):
            raise TypeError("coordinates_u must be VectorQuantity")

    def _check_arg_coordinates_v(self) -> None:
        """Require vector-valued coordinates along the second in-plane axis."""
        if not isinstance(self.coordinates_v, VectorQuantity):
            raise TypeError("coordinates_v must be VectorQuantity")

    def _check_arg_fixed_coordinate(self) -> None:
        """Require a scalar coordinate along the plane-normal axis."""
        if not isinstance(self.fixed_coordinate, ScalarQuantity):
            raise TypeError("fixed_coordinate must be ScalarQuantity")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class Piab3DAnalyticalEigenfunctionPlaneSlice(ResultsObject):
    """Correlate one eigenvector request with a sampled Cartesian plane slice."""

    request: Piab3DAnalyticalEigenfunctionPlaneSliceEvaluationRequest
    coordinates_u: VectorQuantity
    coordinates_v: VectorQuantity
    fixed_coordinate: ScalarQuantity
    amplitudes: MatrixQuantity

    def __post_init__(self) -> None:
        """Check each eigenfunction-plane-slice response argument."""
        self._check_arg_request()
        self._check_arg_coordinates_u()
        self._check_arg_coordinates_v()
        self._check_arg_fixed_coordinate()
        self._check_arg_amplitudes()

    def _check_arg_request(self) -> None:
        """Require the complete eigenfunction plane-slice request."""
        if not isinstance(
            self.request,
            Piab3DAnalyticalEigenfunctionPlaneSliceEvaluationRequest,
        ):
            raise TypeError(
                "request must be "
                "Piab3DAnalyticalEigenfunctionPlaneSliceEvaluationRequest"
            )

    def _check_arg_coordinates_u(self) -> None:
        """Require nonempty model-unit coordinates inside the first plane axis."""
        if not isinstance(self.coordinates_u, VectorQuantity):
            raise TypeError("coordinates_u must be VectorQuantity")
        if self.coordinates_u.unit != self.model.unit_system.length_unit:
            raise ValueError("coordinates_u must use the model length unit")
        if self.coordinates_u.magnitude.size == 0:
            raise ValueError("coordinates_u must be nonempty")
        length_u, _ = self.in_plane_length_magnitudes
        if np.any(self.coordinates_u.magnitude < 0.0) or np.any(
            self.coordinates_u.magnitude > length_u
        ):
            raise ValueError("coordinates_u must lie in its closed box interval")

    def _check_arg_coordinates_v(self) -> None:
        """Require nonempty model-unit coordinates inside the second plane axis."""
        if not isinstance(self.coordinates_v, VectorQuantity):
            raise TypeError("coordinates_v must be VectorQuantity")
        if self.coordinates_v.unit != self.model.unit_system.length_unit:
            raise ValueError("coordinates_v must use the model length unit")
        if self.coordinates_v.magnitude.size == 0:
            raise ValueError("coordinates_v must be nonempty")
        _, length_v = self.in_plane_length_magnitudes
        if np.any(self.coordinates_v.magnitude < 0.0) or np.any(
            self.coordinates_v.magnitude > length_v
        ):
            raise ValueError("coordinates_v must lie in its closed box interval")

    def _check_arg_fixed_coordinate(self) -> None:
        """Require a model-unit coordinate inside the normal-axis interval."""
        if not isinstance(self.fixed_coordinate, ScalarQuantity):
            raise TypeError("fixed_coordinate must be ScalarQuantity")
        if self.fixed_coordinate.unit != self.model.unit_system.length_unit:
            raise ValueError("fixed_coordinate must use the model length unit")
        if not 0.0 <= self.fixed_coordinate.magnitude <= self.normal_length_magnitude:
            raise ValueError("fixed_coordinate must lie in its closed box interval")

    def _check_arg_amplitudes(self) -> None:
        """Require normalized eigenvector amplitudes on the rectangular plane."""
        if not isinstance(self.amplitudes, MatrixQuantity):
            raise TypeError("amplitudes must be MatrixQuantity")
        expected_shape = (
            self.coordinates_u.magnitude.size,
            self.coordinates_v.magnitude.size,
        )
        if self.amplitudes.magnitude.shape != expected_shape:
            raise ValueError("amplitudes must match the rectangular plane")
        length_unit = self.model.unit_system.length_unit
        if isinstance(length_unit, Unitless):
            expected_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(length_unit, PhysicalUnit):
                raise TypeError("physical model requires a physical length unit")
            expected_unit = PhysicalUnit(
                expression=f"({length_unit.expression}) ** -1.5"
            )
        if self.amplitudes.unit != expected_unit:
            raise ValueError("amplitudes must use the normalized eigenfunction unit")

    @property
    def eigenfunction(self) -> Piab3DAnalyticalEigenfunction:
        """Return the stationary eigenvector represented by this plane slice."""
        return self.request.eigenfunction

    @property
    def model(self) -> Piab3D:
        """Return the PIAB3D model that owns the eigenvector."""
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
    def quantum_number_z(self) -> int:
        """Return the eigenvector's z quantum number."""
        return self.eigenfunction.quantum_number_z

    @property
    def plane_normal(self) -> Piab3DPlaneNormal:
        """Return the Cartesian axis fixed by this plane slice."""
        return self.request.plane_normal

    @property
    def in_plane_axis_labels(self) -> tuple[str, str]:
        """Return the Cartesian labels corresponding to ``u`` and ``v``."""
        if self.plane_normal is Piab3DPlaneNormal.X:
            return ("y", "z")
        if self.plane_normal is Piab3DPlaneNormal.Y:
            return ("x", "z")
        return ("x", "y")

    @property
    def in_plane_length_magnitudes(self) -> tuple[float, float]:
        """Return ``u`` and ``v`` lengths in the model length unit."""
        if self.plane_normal is Piab3DPlaneNormal.X:
            return (self.model.length_y, self.model.length_z)
        if self.plane_normal is Piab3DPlaneNormal.Y:
            return (self.model.length_x, self.model.length_z)
        return (self.model.length_x, self.model.length_y)

    @property
    def normal_length_magnitude(self) -> float:
        """Return the fixed-axis length in the model length unit."""
        if self.plane_normal is Piab3DPlaneNormal.X:
            return self.model.length_x
        if self.plane_normal is Piab3DPlaneNormal.Y:
            return self.model.length_y
        return self.model.length_z

    @property
    def normal_quantum_number(self) -> int:
        """Return the quantum number along the fixed-axis direction."""
        if self.plane_normal is Piab3DPlaneNormal.X:
            return self.quantum_number_x
        if self.plane_normal is Piab3DPlaneNormal.Y:
            return self.quantum_number_y
        return self.quantum_number_z

    @property
    def expected_integrated_density(self) -> ScalarQuantity:
        """Return the analytical plane density integrated over its area."""
        normal_length = self.normal_length_magnitude
        phase = (
            self.normal_quantum_number
            * np.pi
            * self.fixed_coordinate.magnitude
            / normal_length
        )
        magnitude = float(2.0 / normal_length * np.sin(phase) ** 2)
        length_unit = self.model.unit_system.length_unit
        if isinstance(length_unit, Unitless):
            density_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(length_unit, PhysicalUnit):
                raise TypeError("physical model requires a physical length unit")
            density_unit = PhysicalUnit(expression=f"({length_unit.expression}) ** -1")
        return ScalarQuantity(magnitude=magnitude, unit=density_unit)

    @property
    def probability_density(self) -> MatrixQuantity:
        r"""Return the sampled probability density on the selected plane."""
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
class Piab3DAnalyticalEigenfunctionPlaneSliceEvaluator:
    """Sample one normalized PIAB3D eigenvector on a Cartesian plane."""

    def action(
        self,
        *,
        request: Piab3DAnalyticalEigenfunctionPlaneSliceEvaluationRequest,
    ) -> Piab3DAnalyticalEigenfunctionPlaneSlice:
        """Return the requested stationary-eigenvector plane slice."""
        if not isinstance(
            request,
            Piab3DAnalyticalEigenfunctionPlaneSliceEvaluationRequest,
        ):
            raise TypeError(
                "request must be "
                "Piab3DAnalyticalEigenfunctionPlaneSliceEvaluationRequest"
            )
        eigenfunction = request.eigenfunction
        model = eigenfunction.model
        length_unit = model.unit_system.length_unit
        selected_u = ModelSystemUnitConverter.convert_vector(
            request.coordinates_u,
            length_unit,
        )
        selected_v = ModelSystemUnitConverter.convert_vector(
            request.coordinates_v,
            length_unit,
        )
        selected_fixed = ModelSystemUnitConverter.convert_scalar(
            request.fixed_coordinate,
            length_unit,
        )
        if selected_u.magnitude.size == 0:
            raise ValueError("coordinates_u must be nonempty")
        if selected_v.magnitude.size == 0:
            raise ValueError("coordinates_v must be nonempty")

        if request.plane_normal is Piab3DPlaneNormal.X:
            length_u, length_v, normal_length = (
                model.length_y,
                model.length_z,
                model.length_x,
            )
            quantum_u, quantum_v, normal_quantum_number = (
                eigenfunction.quantum_number_y,
                eigenfunction.quantum_number_z,
                eigenfunction.quantum_number_x,
            )
        elif request.plane_normal is Piab3DPlaneNormal.Y:
            length_u, length_v, normal_length = (
                model.length_x,
                model.length_z,
                model.length_y,
            )
            quantum_u, quantum_v, normal_quantum_number = (
                eigenfunction.quantum_number_x,
                eigenfunction.quantum_number_z,
                eigenfunction.quantum_number_y,
            )
        else:
            length_u, length_v, normal_length = (
                model.length_x,
                model.length_y,
                model.length_z,
            )
            quantum_u, quantum_v, normal_quantum_number = (
                eigenfunction.quantum_number_x,
                eigenfunction.quantum_number_y,
                eigenfunction.quantum_number_z,
            )
        if np.any(selected_u.magnitude < 0.0) or np.any(
            selected_u.magnitude > length_u
        ):
            raise ValueError("coordinates_u must lie in its closed box interval")
        if np.any(selected_v.magnitude < 0.0) or np.any(
            selected_v.magnitude > length_v
        ):
            raise ValueError("coordinates_v must lie in its closed box interval")
        if not 0.0 <= selected_fixed.magnitude <= normal_length:
            raise ValueError("fixed_coordinate must lie in its closed box interval")

        u_factor = np.sin(quantum_u * np.pi * selected_u.magnitude / length_u)
        v_factor = np.sin(quantum_v * np.pi * selected_v.magnitude / length_v)
        normal_factor = np.sin(
            normal_quantum_number * np.pi * selected_fixed.magnitude / normal_length
        )
        values = (
            np.sqrt(8.0 / (model.length_x * model.length_y * model.length_z))
            * normal_factor
            * np.outer(u_factor, v_factor)
        )
        if isinstance(length_unit, Unitless):
            amplitude_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(length_unit, PhysicalUnit):
                raise TypeError("physical model requires a physical length unit")
            amplitude_unit = PhysicalUnit(
                expression=f"({length_unit.expression}) ** -1.5"
            )
        return Piab3DAnalyticalEigenfunctionPlaneSlice(
            request=request,
            coordinates_u=selected_u,
            coordinates_v=selected_v,
            fixed_coordinate=selected_fixed,
            amplitudes=MatrixQuantity(
                magnitude=values,
                unit=amplitude_unit,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class Piab3DAnalyticalEigenfunction(
    ActionizedDataObject[
        Piab3DAnalyticalEigenfunctionPlaneSliceEvaluationRequest,
        Piab3DAnalyticalEigenfunctionPlaneSlice,
    ]
):
    """Represent one stationary PIAB3D Hamiltonian eigenvector."""

    model: Piab3D
    quantum_number_x: int
    quantum_number_y: int
    quantum_number_z: int
    actionizer: Piab3DAnalyticalEigenfunctionPlaneSliceEvaluator = field(
        default_factory=Piab3DAnalyticalEigenfunctionPlaneSliceEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check each analytical-eigenfunction argument."""
        self._check_arg_model()
        self._check_arg_quantum_number_x()
        self._check_arg_quantum_number_y()
        self._check_arg_quantum_number_z()

    def _check_arg_model(self) -> None:
        """Require the cuboid PIAB3D Hamiltonian model."""
        if not isinstance(self.model, Piab3D):
            raise TypeError("model must be Piab3D")

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

    def _check_arg_quantum_number_z(self) -> None:
        """Require a positive built-in integer z quantum number."""
        if type(self.quantum_number_z) is not int:
            raise TypeError("quantum_number_z must be a built-in integer")
        if self.quantum_number_z <= 0:
            raise ValueError("quantum_number_z must be positive")

    def evaluate_plane_slice(
        self,
        *,
        plane_normal: Piab3DPlaneNormal,
        coordinates_u: VectorQuantity,
        coordinates_v: VectorQuantity,
        fixed_coordinate: ScalarQuantity,
    ) -> Piab3DAnalyticalEigenfunctionPlaneSlice:
        """Sample this eigenvector on one Cartesian plane."""
        return self._respond(
            request=Piab3DAnalyticalEigenfunctionPlaneSliceEvaluationRequest(
                eigenfunction=self,
                plane_normal=plane_normal,
                coordinates_u=coordinates_u,
                coordinates_v=coordinates_v,
                fixed_coordinate=fixed_coordinate,
            )
        )
