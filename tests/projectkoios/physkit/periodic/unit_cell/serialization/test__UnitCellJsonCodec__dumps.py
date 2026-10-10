from __future__ import annotations

import json

import pytest

from projectkoios.physkit.periodic.lattice.lattice3d import DirectLattice3D
from projectkoios.physkit.periodic.unit_cell.base import (
    AtomicBasis,
    ConventionalUnitCell,
    PrimitiveUnitCell,
    UnitCell,
)
from projectkoios.physkit.periodic.unit_cell.serialization import UnitCellJsonCodec
from projectkoios.physkit.units.quantities import PhysicalUnit, ScalarQuantity


def test_emits_deterministic_version_one_payload(
    fcc_direct_lattice: DirectLattice3D,
    silicon_basis: AtomicBasis,
) -> None:
    cell = PrimitiveUnitCell(
        direct_lattice=fcc_direct_lattice,
        lattice_parameter=ScalarQuantity(5.43, PhysicalUnit("angstrom")),
        atomic_basis=silicon_basis,
    )

    first = UnitCellJsonCodec().dumps(
        cell,
        structure_id="Si.PrimitiveUnitCell",
    )
    second = UnitCellJsonCodec().dumps(
        cell,
        structure_id="Si.PrimitiveUnitCell",
    )
    payload = json.loads(first)

    assert first == second
    assert first.endswith("\n")
    assert payload["schema_version"] == 1
    assert payload["structure_id"] == "Si.PrimitiveUnitCell"
    assert payload["representation"] == "primitive"
    assert payload["lattice"]["A"]["vector_axis"] == "columns"


@pytest.mark.parametrize(
    ("cell_type", "representation", "structure_id"),
    [
        (PrimitiveUnitCell, "primitive", "Si.PrimitiveUnitCell"),
        (ConventionalUnitCell, "conventional", "Si.ConventionalUnitCell"),
    ],
)
def test__UnitCellJsonCodec_dumps__existing_representations__preserve_bytes(
    cell_type: type[UnitCell],
    representation: str,
    structure_id: str,
    fcc_direct_lattice: DirectLattice3D,
    silicon_basis: AtomicBasis,
) -> None:
    cell = cell_type(
        direct_lattice=fcc_direct_lattice,
        lattice_parameter=ScalarQuantity(5.43, PhysicalUnit("angstrom")),
        atomic_basis=silicon_basis,
    )
    expected = f'''{{
  "atoms": [
    {{
      "fractional_position": [
        0.0,
        0.0,
        0.0
      ],
      "symbol": "Si"
    }},
    {{
      "fractional_position": [
        0.25,
        0.25,
        0.25
      ],
      "symbol": "Si"
    }}
  ],
  "lattice": {{
    "A": {{
      "a1": [
        0.5,
        0.5,
        0.0
      ],
      "a2": [
        0.5,
        0.0,
        0.5
      ],
      "a3": [
        0.0,
        0.5,
        0.5
      ],
      "vector_axis": "columns"
    }},
    "lattice_parameter": {{
      "magnitude": 5.43,
      "unit": "angstrom"
    }}
  }},
  "representation": "{representation}",
  "schema_version": 1,
  "structure_id": "{structure_id}"
}}
'''

    assert UnitCellJsonCodec().dumps(cell, structure_id=structure_id) == expected


def test__UnitCellJsonCodec_dumps__base_unit_cell__emits_deterministic_schema(
    fcc_direct_lattice: DirectLattice3D,
    silicon_basis: AtomicBasis,
) -> None:
    cell = UnitCell(
        direct_lattice=fcc_direct_lattice,
        lattice_parameter=ScalarQuantity(5.43, PhysicalUnit("angstrom")),
        atomic_basis=silicon_basis,
    )
    codec = UnitCellJsonCodec()

    first = codec.dumps(cell, structure_id="Si.UnitCell")
    second = codec.dumps(cell, structure_id="Si.UnitCell")
    payload = json.loads(first)

    assert first == second
    assert first.endswith("\n")
    assert payload["representation"] == "unit-cell"
    assert payload["schema_version"] == 1
    assert payload["lattice"]["A"]["vector_axis"] == "columns"
    assert [atom["symbol"] for atom in payload["atoms"]] == ["Si", "Si"]
