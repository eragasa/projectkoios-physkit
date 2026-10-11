"""Deterministic encoding tests for embedded 2D lattice JSON."""

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.embedding import EmbeddedDirectLattice2D
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


def test__dumps__emits_deterministic_physical_schema_golden_text() -> None:
    codec = EmbeddedDirectLattice2DJsonCodec()
    embedded = _embedded(PhysicalUnit("angstrom"))
    expected = (
        '{"ambient_basis":{"a1":[2.0,0.3,1.1],'
        '"a2":[0.4,1.7,-0.2],'
        '"unit":{"expression":"angstrom","kind":"physical"},'
        '"vector_axis":"columns"},'
        '"representation":"embedded-direct-lattice-2d",'
        '"schema_version":1}\n'
    )

    assert codec.dumps(embedded) == expected
    assert codec.dumps(embedded) == codec.dumps(embedded)


def test__dumps__emits_explicit_unitless_discriminator() -> None:
    content = EmbeddedDirectLattice2DJsonCodec().dumps(_embedded(Unitless()))

    assert '"unit":{"kind":"unitless"}' in content
    assert '"expression"' not in content
    assert content.endswith("\n")


def test__dumps__rejects_nonexact_nominal_type() -> None:
    class SpecializedEmbeddedDirectLattice2D(EmbeddedDirectLattice2D):
        pass

    specialized = SpecializedEmbeddedDirectLattice2D(
        ambient_basis=_embedded(Unitless()).ambient_basis
    )

    with pytest.raises(
        EmbeddedDirectLattice2DSerializationError,
        match="exactly EmbeddedDirectLattice2D",
    ):
        EmbeddedDirectLattice2DJsonCodec().dumps(specialized)
