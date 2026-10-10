"""Versioned deterministic JSON serialization for periodic unit cells."""

from __future__ import annotations

import json

import numpy as np

from projectkoios.physkit.periodic.lattice.lattice3d import DirectLattice3D
from projectkoios.physkit.periodic.unit_cell.base import (
    Atom,
    AtomicBasis,
    ConventionalUnitCell,
    PrimitiveUnitCell,
    UnitCell,
)
from projectkoios.physkit.periodic.unit_cell.errors import UnitCellSerializationError
from projectkoios.physkit.units.quantities import (
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

_SCHEMA_VERSION = 1
_VECTOR_AXIS = "columns"
_REPRESENTATION_BY_TYPE = {
    ConventionalUnitCell: "conventional",
    PrimitiveUnitCell: "primitive",
    UnitCell: "unit-cell",
}
_CELL_TYPE_BY_REPRESENTATION: dict[str, type[UnitCell]] = {
    value: key for key, value in _REPRESENTATION_BY_TYPE.items()
}


class UnitCellJsonCodec:
    """Encode and decode the version-one deterministic unit-cell JSON schema."""

    __slots__ = ()

    def dumps(self, unit_cell: UnitCell, *, structure_id: str) -> str:
        """Return deterministic JSON for one exact supported unit-cell type."""
        if type(structure_id) is not str or not structure_id:
            raise UnitCellSerializationError("structure_id must be a nonempty string")
        representation = _REPRESENTATION_BY_TYPE.get(type(unit_cell))
        if representation is None:
            raise UnitCellSerializationError(
                "unit_cell must be exactly UnitCell, ConventionalUnitCell, "
                "or PrimitiveUnitCell"
            )
        length_unit = unit_cell.lattice_parameter.unit
        if type(length_unit) is not PhysicalUnit:
            raise UnitCellSerializationError(
                "unit_cell.lattice_parameter must have a physical unit"
            )
        payload = {
            "atoms": [
                {
                    "fractional_position": atom.position_fractional.magnitude.tolist(),
                    "symbol": atom.symbol,
                }
                for atom in unit_cell.atomic_basis.atoms
            ],
            "lattice": {
                "A": {
                    "a1": unit_cell.a1.magnitude.tolist(),
                    "a2": unit_cell.a2.magnitude.tolist(),
                    "a3": unit_cell.a3.magnitude.tolist(),
                    "vector_axis": _VECTOR_AXIS,
                },
                "lattice_parameter": {
                    "magnitude": unit_cell.lattice_parameter.magnitude,
                    "unit": length_unit.expression,
                },
            },
            "representation": representation,
            "schema_version": _SCHEMA_VERSION,
            "structure_id": structure_id,
        }
        return json.dumps(payload, indent=2, sort_keys=True) + "\n"

    def loads(self, content: str, *, expected_structure_id: str) -> UnitCell:
        """Return one validated unit cell of the explicitly encoded exact type."""
        if type(content) is not str:
            raise UnitCellSerializationError("content must be a string")
        if type(expected_structure_id) is not str or not expected_structure_id:
            raise UnitCellSerializationError(
                "expected_structure_id must be a nonempty string"
            )
        try:
            payload = json.loads(content)
        except json.JSONDecodeError as error:
            raise UnitCellSerializationError(
                "unit-cell content is not valid JSON"
            ) from error
        if not isinstance(payload, dict):
            raise UnitCellSerializationError("unit-cell document must be an object")
        if payload.get("schema_version") != _SCHEMA_VERSION:
            raise UnitCellSerializationError("unsupported unit-cell JSON schema")
        if payload.get("structure_id") != expected_structure_id:
            raise UnitCellSerializationError("structure record identity does not match")

        representation = _require_string(payload, "representation")
        cell_type = _CELL_TYPE_BY_REPRESENTATION.get(representation)
        if cell_type is None:
            raise UnitCellSerializationError(
                "representation must be 'conventional', 'primitive', or 'unit-cell'"
            )
        lattice = _require_mapping(payload, "lattice")
        a_matrix = _require_mapping(lattice, "A", path="lattice")
        if _require_string(a_matrix, "vector_axis", path="lattice.A") != _VECTOR_AXIS:
            raise UnitCellSerializationError(
                "lattice.A.vector_axis must declare vectors as columns"
            )
        raw_atoms = payload.get("atoms")
        if not isinstance(raw_atoms, list) or not raw_atoms:
            raise UnitCellSerializationError("atoms must be a nonempty array")
        lattice_parameter = _require_mapping(
            lattice,
            "lattice_parameter",
            path="lattice",
        )

        try:
            return cell_type(
                direct_lattice=DirectLattice3D(
                    a1=np.array(_decode_triplet(a_matrix, "a1", path="lattice.A")),
                    a2=np.array(_decode_triplet(a_matrix, "a2", path="lattice.A")),
                    a3=np.array(_decode_triplet(a_matrix, "a3", path="lattice.A")),
                ),
                lattice_parameter=ScalarQuantity(
                    magnitude=_require_number(
                        lattice_parameter.get("magnitude"),
                        "lattice.lattice_parameter.magnitude",
                    ),
                    unit=PhysicalUnit(
                        _require_string(
                            lattice_parameter,
                            "unit",
                            path="lattice.lattice_parameter",
                        )
                    ),
                ),
                atomic_basis=AtomicBasis(
                    atoms=tuple(
                        _decode_atom(value, index=index)
                        for index, value in enumerate(raw_atoms)
                    )
                ),
            )
        except UnitCellSerializationError:
            raise
        except (TypeError, ValueError) as error:
            raise UnitCellSerializationError(
                "decoded unit-cell values do not form a valid unit cell"
            ) from error


def _decode_atom(value: object, *, index: int) -> Atom:
    path = f"atoms[{index}]"
    if not isinstance(value, dict):
        raise UnitCellSerializationError(f"{path} must be an object")
    return Atom(
        symbol=_require_string(value, "symbol", path=path),
        position_fractional=VectorQuantity(
            magnitude=np.array(
                _decode_triplet(value, "fractional_position", path=path)
            ),
            unit=Unitless(),
        ),
    )


def _decode_triplet(
    mapping: dict[str, object],
    key: str,
    *,
    path: str,
) -> tuple[float, float, float]:
    value = mapping.get(key)
    qualified = f"{path}.{key}"
    if not isinstance(value, list) or len(value) != 3:
        raise UnitCellSerializationError(
            f"{qualified} must contain three finite numbers"
        )
    numbers = tuple(_require_number(item, qualified) for item in value)
    return numbers[0], numbers[1], numbers[2]


def _require_mapping(
    mapping: dict[str, object],
    key: str,
    *,
    path: str = "",
) -> dict[str, object]:
    value = mapping.get(key)
    qualified = f"{path}.{key}" if path else key
    if not isinstance(value, dict):
        raise UnitCellSerializationError(f"{qualified} must be an object")
    return value


def _require_string(
    mapping: dict[str, object],
    key: str,
    *,
    path: str = "",
) -> str:
    value = mapping.get(key)
    qualified = f"{path}.{key}" if path else key
    if type(value) is not str or not value or value != value.strip():
        raise UnitCellSerializationError(
            f"{qualified} must be a nonempty stripped string"
        )
    return value


def _require_number(value: object, path: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise UnitCellSerializationError(f"{path} must be a number")
    number = float(value)
    if not np.isfinite(number):
        raise UnitCellSerializationError(f"{path} must be finite")
    return number
