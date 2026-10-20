"""Unit and identity correlation across finite-periodic matrix operations."""

import numpy as np
import pytest

from projectkoios.numerics.finite_periodic.meshes import CenteredUniformMesh
from projectkoios.physkit.solidstate.tb.finite_periodic import (
    FinitePeriodicMatrixFourierEvaluation,
    FinitePeriodicMatrixFourierSeries,
    FinitePeriodicMatrixFourierTensorEvaluator,
    FinitePeriodicMatrixFourierTransformer,
    FinitePeriodicMatrixFourierTransformResult,
    FinitePeriodicMatrixHermiticityAnalyzer,
    FinitePeriodicMatrixHermiticityDiagnostics,
    FinitePeriodicMatrixParsevalAnalyzer,
    FinitePeriodicMatrixParsevalDiagnostics,
    FinitePeriodicMatrixSamples,
    diagnostics,
)
from projectkoios.physkit.units import PhysicalUnit, Unitless


def _hermitian_samples(mesh: CenteredUniformMesh) -> np.ndarray:
    axis = np.asarray(
        [-0.5 + (index + mesh.offsets[0]) / mesh.size for index in range(mesh.size)],
        dtype=np.float64,
    )
    values = np.empty((mesh.size, 2, 2), dtype=np.complex128)
    values[:, 0, 0] = 1.0 + np.cos(2.0 * np.pi * axis)
    values[:, 1, 1] = 2.0 - np.cos(2.0 * np.pi * axis)
    values[:, 0, 1] = 0.25j * np.sin(2.0 * np.pi * axis)
    values[:, 1, 0] = values[:, 0, 1].conjugate()
    return values


@pytest.mark.parametrize("size", [5, 6])
def test_transform_evaluation_and_diagnostics_preserve_energy_unit(size: int) -> None:
    mesh = CenteredUniformMesh(1, size, (0.0,))
    samples = _hermitian_samples(mesh)
    unit = PhysicalUnit("electron_volt")
    source = FinitePeriodicMatrixSamples.from_array(
        mesh=mesh,
        matrix_samples=samples,
        energy_unit=unit,
        energy_reference="declared-zero",
        basis_identifier="ordered-two-state-basis",
    )

    transform = FinitePeriodicMatrixFourierTransformer().transform(source)
    axis = np.asarray(
        [-0.5 + index / size for index in range(size)],
        dtype=np.float64,
    )
    evaluation = FinitePeriodicMatrixFourierTensorEvaluator().evaluate(
        transform.series,
        (axis,),
    )
    parseval = FinitePeriodicMatrixParsevalAnalyzer().analyze(transform)
    hermiticity = FinitePeriodicMatrixHermiticityAnalyzer().analyze(transform)

    expected_modes = np.arange(-(size // 2), (size - 1) // 2 + 1)
    assert transform.source is source
    assert transform.series.source_mesh is mesh
    assert np.array_equal(transform.series.integer_modes[:, 0], expected_modes)
    zero_index = int(np.flatnonzero(expected_modes == 0)[0])
    assert np.allclose(
        transform.series.matrix_coefficients[zero_index],
        np.mean(samples, axis=0),
        atol=2.0e-15,
        rtol=0.0,
    )
    assert transform.backend_identifier == "numpy.fft.fftn"
    assert transform.backend_version
    assert np.allclose(evaluation.matrix_values, samples, atol=2.0e-15, rtol=0.0)
    assert evaluation.backend_identifier == "numpy.einsum.tensor_product"
    assert evaluation.backend_version
    assert transform.series.energy_unit is unit
    assert evaluation.energy_unit is unit
    assert evaluation.energy_reference == "declared-zero"
    assert evaluation.basis_identifier == "ordered-two-state-basis"
    squared = PhysicalUnit("(electron_volt) ** 2")
    assert parseval.sample_mean_squared_frobenius_norm.unit == squared
    assert parseval.coefficient_sum_squared_frobenius_norm.unit == squared
    assert parseval.absolute_defect.unit == squared
    assert hermiticity.twisted_cyclic_maximum_frobenius_defect.unit is unit
    assert hermiticity.ordinary_maximum_frobenius_defect.unit is unit
    assert hermiticity.has_nyquist_boundary is (size % 2 == 0)
    assert hermiticity.ordinary_comparison_covers_complete_inventory is (size % 2 == 1)


def test_computed_results_require_their_owning_actions() -> None:
    mesh = CenteredUniformMesh(1, 5, (0.0,))
    source = FinitePeriodicMatrixSamples.from_array(
        mesh=mesh,
        matrix_samples=_hermitian_samples(mesh),
        energy_unit=PhysicalUnit("eV"),
        energy_reference="zero",
        basis_identifier="basis",
    )
    transform = FinitePeriodicMatrixFourierTransformer().transform(source)
    evaluation = FinitePeriodicMatrixFourierTensorEvaluator().evaluate(
        transform.series,
        (np.asarray([0.0], dtype=np.float64),),
    )
    parseval = FinitePeriodicMatrixParsevalAnalyzer().analyze(transform)
    hermiticity = FinitePeriodicMatrixHermiticityAnalyzer().analyze(transform)

    with pytest.raises(TypeError, match="owning action"):
        FinitePeriodicMatrixFourierTransformResult(
            source=transform.source,
            series=transform.series,
            backend_identifier=transform.backend_identifier,
            backend_version=transform.backend_version,
        )
    with pytest.raises(TypeError, match="owning action"):
        FinitePeriodicMatrixFourierEvaluation(
            series=evaluation.series,
            coordinate_axes=evaluation.coordinate_axes,
            matrix_values=evaluation.matrix_values,
            maximum_working_elements=evaluation.maximum_working_elements,
            backend_identifier=evaluation.backend_identifier,
            backend_version=evaluation.backend_version,
        )
    with pytest.raises(TypeError, match="owning action"):
        FinitePeriodicMatrixParsevalDiagnostics(
            transform=parseval.transform,
            sample_mean_squared_frobenius_norm=(
                parseval.sample_mean_squared_frobenius_norm
            ),
            coefficient_sum_squared_frobenius_norm=(
                parseval.coefficient_sum_squared_frobenius_norm
            ),
            absolute_defect=parseval.absolute_defect,
        )
    with pytest.raises(TypeError, match="owning action"):
        FinitePeriodicMatrixHermiticityDiagnostics(
            transform=hermiticity.transform,
            twisted_cyclic_maximum_frobenius_defect=(
                hermiticity.twisted_cyclic_maximum_frobenius_defect
            ),
            ordinary_maximum_frobenius_defect=(
                hermiticity.ordinary_maximum_frobenius_defect
            ),
        )


def test_each_diagnostic_owner_runs_once(monkeypatch: pytest.MonkeyPatch) -> None:
    mesh = CenteredUniformMesh(1, 5, (0.0,))
    source = FinitePeriodicMatrixSamples.from_array(
        mesh=mesh,
        matrix_samples=_hermitian_samples(mesh),
        energy_unit=PhysicalUnit("eV"),
        energy_reference="zero",
        basis_identifier="basis",
    )
    transform = FinitePeriodicMatrixFourierTransformer().transform(source)
    calls = {"parseval": 0, "hermiticity": 0}
    parseval_owner = diagnostics.NumpyFinitePeriodicParsevalAnalyzer
    hermiticity_owner = diagnostics.NumpyCenteredUniformMeshHermiticityAnalyzer

    class CountingParsevalAnalyzer:
        def analyze(self, *args: object) -> object:
            calls["parseval"] += 1
            return parseval_owner().analyze(*args)

    class CountingHermiticityAnalyzer:
        def analyze(self, *args: object) -> object:
            calls["hermiticity"] += 1
            return hermiticity_owner().analyze(*args)

    monkeypatch.setattr(
        diagnostics,
        "NumpyFinitePeriodicParsevalAnalyzer",
        CountingParsevalAnalyzer,
    )
    monkeypatch.setattr(
        diagnostics,
        "NumpyCenteredUniformMeshHermiticityAnalyzer",
        CountingHermiticityAnalyzer,
    )

    FinitePeriodicMatrixParsevalAnalyzer().analyze(transform)
    FinitePeriodicMatrixHermiticityAnalyzer().analyze(transform)

    assert calls == {"parseval": 1, "hermiticity": 1}


def test_selected_series_retains_units_and_explicit_inventory_order() -> None:
    mesh = CenteredUniformMesh(1, 5, (0.0,))
    source = FinitePeriodicMatrixSamples.from_array(
        mesh=mesh,
        matrix_samples=_hermitian_samples(mesh),
        energy_unit=PhysicalUnit("eV"),
        energy_reference="zero",
        basis_identifier="basis",
    )
    complete = FinitePeriodicMatrixFourierTransformer().transform(source).series

    selected = complete.take((1, 2, 3))

    assert np.array_equal(selected.integer_modes[:, 0], np.asarray([-1, 0, 1]))
    assert selected.energy_unit is complete.energy_unit
    assert selected.energy_reference == complete.energy_reference
    assert selected.basis_identifier == complete.basis_identifier
    for invalid in ((2, 1), (1, 1), (-1,), (complete.term_count,)):
        with pytest.raises(ValueError):
            complete.take(invalid)
    with pytest.raises(TypeError, match="tuple"):
        complete.take([1, 2])  # type: ignore[arg-type]


def test_empty_selected_series_evaluates_as_exact_zero_with_unit() -> None:
    mesh = CenteredUniformMesh(1, 5, (0.0,))
    source = FinitePeriodicMatrixSamples.from_array(
        mesh=mesh,
        matrix_samples=_hermitian_samples(mesh),
        energy_unit=Unitless(),
        energy_reference="dimensionless-zero",
        basis_identifier="basis",
    )
    empty = FinitePeriodicMatrixFourierTransformer().transform(source).series.take(())
    axis = np.asarray([-0.25, 0.125], dtype=np.float64)

    evaluation = FinitePeriodicMatrixFourierTensorEvaluator().evaluate(empty, (axis,))

    assert evaluation.matrix_values.shape == (2, 2, 2)
    assert np.array_equal(evaluation.matrix_values, np.zeros((2, 2, 2)))
    assert isinstance(evaluation.energy_unit, Unitless)


def test_evaluation_copies_and_freezes_explicit_target_axes() -> None:
    mesh = CenteredUniformMesh(1, 5, (0.0,))
    source = FinitePeriodicMatrixSamples.from_array(
        mesh=mesh,
        matrix_samples=_hermitian_samples(mesh),
        energy_unit=PhysicalUnit("eV"),
        energy_reference="zero",
        basis_identifier="basis",
    )
    series = FinitePeriodicMatrixFourierTransformer().transform(source).series
    axis = np.asarray([-0.25, 0.125], dtype=np.float64)

    evaluation = FinitePeriodicMatrixFourierTensorEvaluator().evaluate(series, (axis,))
    axis[0] = 99.0

    assert np.array_equal(evaluation.coordinate_axes[0], [-0.25, 0.125])
    assert not evaluation.coordinate_axes[0].flags.writeable
    assert not evaluation.matrix_values.flags.writeable


@pytest.mark.parametrize(
    "modes, message",
    [
        (np.asarray([[3]], dtype=np.int64), "source mesh inventory"),
        (np.asarray([[0], [0]], dtype=np.int64), "ordered and unique"),
        (np.asarray([[1], [-1]], dtype=np.int64), "ordered and unique"),
    ],
)
def test_series_rejects_invalid_physical_inventories(
    modes: np.ndarray,
    message: str,
) -> None:
    mesh = CenteredUniformMesh(1, 5, (0.0,))
    coefficients = np.ones((modes.shape[0], 1, 1), dtype=np.complex128)

    with pytest.raises(ValueError, match=message):
        FinitePeriodicMatrixFourierSeries.from_arrays(
            source_mesh=mesh,
            integer_modes=modes,
            matrix_coefficients=coefficients,
            energy_unit=PhysicalUnit("eV"),
            energy_reference="zero",
            basis_identifier="basis",
        )


@pytest.mark.parametrize(
    "values, exception",
    [
        (np.zeros((5, 1, 1), dtype=np.float64), TypeError),
        (
            np.asarray([[[np.nan + 0.0j]]] * 5, dtype=np.complex128),
            ValueError,
        ),
    ],
)
def test_source_rejects_wrong_dtype_and_nonfinite_arrays(
    values: np.ndarray,
    exception: type[Exception],
) -> None:
    with pytest.raises(exception):
        FinitePeriodicMatrixSamples.from_array(
            mesh=CenteredUniformMesh(1, 5, (0.0,)),
            matrix_samples=values,
            energy_unit=PhysicalUnit("eV"),
            energy_reference="zero",
            basis_identifier="basis",
        )


def test_source_factory_copies_and_freezes_samples() -> None:
    mesh = CenteredUniformMesh(1, 5, (0.0,))
    values = _hermitian_samples(mesh)
    source = FinitePeriodicMatrixSamples.from_array(
        mesh=mesh,
        matrix_samples=values,
        energy_unit=PhysicalUnit("eV"),
        energy_reference="zero",
        basis_identifier="basis",
    )

    values[0, 0, 0] = 99.0

    assert source.matrix_samples[0, 0, 0] != 99.0
    assert not source.matrix_samples.flags.writeable
