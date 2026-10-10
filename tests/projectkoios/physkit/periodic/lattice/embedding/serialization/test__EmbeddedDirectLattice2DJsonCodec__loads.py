"""Strict decoding tests for embedded 2D lattice JSON."""

from __future__ import annotations

import json

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.embedding import (
    EmbeddedDirectLattice2D,
    EmbeddedDirectLattice2DAdaptationRequest,
    EmbeddedDirectLattice2DAdapter,
)
from projectkoios.physkit.periodic.lattice.embedding.serialization import (
    EmbeddedDirectLattice2DJsonCodec,
    EmbeddedDirectLattice2DSerializationError,
)
from projectkoios.physkit.units import MatrixQuantity, PhysicalUnit, Unitless


def _embedded(unit: PhysicalUnit | Unitless) -> EmbeddedDirectLattice2D:
    return EmbeddedDirectLattice2D(
        ambient_basis=MatrixQuantity(
            np.array(
                [
                    [2.0, 0.4],
                    [0.3, 1.7],
                    [1.1, -0.2],
                ]
            ),
            unit,
        )
    )


def _payload(
    unit: PhysicalUnit | Unitless | None = None,
) -> dict[str, object]:
    selected_unit = Unitless() if unit is None else unit
    return json.loads(
        EmbeddedDirectLattice2DJsonCodec().dumps(_embedded(selected_unit))
    )


def test__loads__round_trip_preserves_exact_type_values_units_and_bytes() -> None:
    codec = EmbeddedDirectLattice2DJsonCodec()
    for source in (
        _embedded(Unitless()),
        _embedded(PhysicalUnit("angstrom")),
    ):
        content = codec.dumps(source)
        restored = codec.loads(content)

        assert type(restored) is EmbeddedDirectLattice2D
        assert restored.ambient_basis.unit == source.ambient_basis.unit
        np.testing.assert_array_equal(
            restored.ambient_basis.magnitude,
            source.ambient_basis.magnitude,
        )
        assert codec.dumps(restored) == content


def test__loads__reconstructs_derived_geometry_from_defining_basis() -> None:
    codec = EmbeddedDirectLattice2DJsonCodec()
    source = _embedded(PhysicalUnit("angstrom"))
    restored = codec.loads(codec.dumps(source))
    source_geometry = EmbeddedDirectLattice2DAdapter().action(
        request=EmbeddedDirectLattice2DAdaptationRequest(embedded_lattice=source)
    )
    restored_geometry = EmbeddedDirectLattice2DAdapter().action(
        request=EmbeddedDirectLattice2DAdaptationRequest(embedded_lattice=restored)
    )

    np.testing.assert_array_equal(
        restored_geometry.induced_metric.magnitude,
        source_geometry.induced_metric.magnitude,
    )
    np.testing.assert_array_equal(
        restored_geometry.inverse_induced_metric.magnitude,
        source_geometry.inverse_induced_metric.magnitude,
    )
    np.testing.assert_array_equal(
        restored_geometry.dual_reciprocal_basis.magnitude,
        source_geometry.dual_reciprocal_basis.magnitude,
    )
    np.testing.assert_array_equal(
        restored_geometry.ambient_projector.magnitude,
        source_geometry.ambient_projector.magnitude,
    )
    np.testing.assert_array_equal(
        restored_geometry.orthonormal_frame.magnitude,
        source_geometry.orthonormal_frame.magnitude,
    )


@pytest.mark.parametrize(
    ("path", "key"),
    [
        ((), "representation"),
        (("ambient_basis",), "a2"),
        (("ambient_basis", "unit"), "kind"),
    ],
)
def test__loads__rejects_missing_members(
    path: tuple[str, ...],
    key: str,
) -> None:
    payload = _payload()
    target = payload
    for component in path:
        target = target[component]  # type: ignore[assignment,index]
    del target[key]  # type: ignore[index]

    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="missing",
    ):
        EmbeddedDirectLattice2DJsonCodec().loads(json.dumps(payload))


@pytest.mark.parametrize(
    "path",
    [(), ("ambient_basis",), ("ambient_basis", "unit")],
)
def test__loads__rejects_unknown_members(path: tuple[str, ...]) -> None:
    payload = _payload()
    target = payload
    for component in path:
        target = target[component]  # type: ignore[assignment,index]
    target["derived_metric"] = [[1.0, 0.0], [0.0, 1.0]]  # type: ignore[index]

    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="unknown",
    ):
        EmbeddedDirectLattice2DJsonCodec().loads(json.dumps(payload))


@pytest.mark.parametrize(
    "content",
    [
        '{"schema_version":1,"schema_version":1}',
        (
            '{"ambient_basis":{"a1":[1,0,0],"a1":[1,0,0],'
            '"a2":[0,1,0],"unit":{"kind":"unitless"},'
            '"vector_axis":"columns"},'
            '"representation":"embedded-direct-lattice-2d",'
            '"schema_version":1}'
        ),
        (
            '{"ambient_basis":{"a1":[1,0,0],"a2":[0,1,0],'
            '"unit":{"kind":"unitless","kind":"unitless"},'
            '"vector_axis":"columns"},'
            '"representation":"embedded-direct-lattice-2d",'
            '"schema_version":1}'
        ),
    ],
)
def test__loads__rejects_duplicate_members(content: str) -> None:
    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="duplicate JSON object key",
    ):
        EmbeddedDirectLattice2DJsonCodec().loads(content)


@pytest.mark.parametrize("token", ["NaN", "Infinity", "-Infinity", "1e999"])
def test__loads__rejects_nonfinite_numbers(token: str) -> None:
    content = EmbeddedDirectLattice2DJsonCodec().dumps(_embedded(Unitless()))
    content = content.replace("2.0", token, 1)

    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="finite|constant|strict valid JSON",
    ):
        EmbeddedDirectLattice2DJsonCodec().loads(content)


@pytest.mark.parametrize("invalid", [True, "2.0", None, {}, []])
def test__loads__rejects_non_numeric_vector_values(invalid: object) -> None:
    payload = _payload()
    payload["ambient_basis"]["a1"][0] = invalid  # type: ignore[index]

    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="JSON number",
    ):
        EmbeddedDirectLattice2DJsonCodec().loads(json.dumps(payload))


def test__loads__rejects_schema_representation_axis_and_unit_variants() -> None:
    codec = EmbeddedDirectLattice2DJsonCodec()

    payload = _payload()
    payload["schema_version"] = True
    with pytest.raises(EmbeddedDirectLattice2DSerializationError, match="integer"):
        codec.loads(json.dumps(payload))

    payload = _payload()
    payload["representation"] = "embedded-plane"
    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="representation",
    ):
        codec.loads(json.dumps(payload))

    payload = _payload()
    payload["ambient_basis"]["vector_axis"] = "rows"  # type: ignore[index]
    with pytest.raises(EmbeddedDirectLattice2DSerializationError, match="columns"):
        codec.loads(json.dumps(payload))

    payload = _payload()
    payload["ambient_basis"]["unit"] = {  # type: ignore[index]
        "kind": "physical",
        "expression": "electron_volt",
    }
    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="valid embedded lattice",
    ):
        codec.loads(json.dumps(payload))


def test__loads__rejects_invalid_rank_malformed_json_and_resource_excess() -> None:
    codec = EmbeddedDirectLattice2DJsonCodec()
    payload = _payload()
    payload["ambient_basis"]["a2"] = [4.0, 0.6, 2.2]  # type: ignore[index]
    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="valid embedded lattice",
    ):
        codec.loads(json.dumps(payload))

    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="strict valid JSON",
    ):
        codec.loads("{not-json")

    oversized = " " * (codec.MAX_CONTENT_CHARACTERS + 1)
    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="resource bound",
    ):
        codec.loads(oversized)
