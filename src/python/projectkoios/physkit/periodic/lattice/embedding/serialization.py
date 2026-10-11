r"""Strict versioned JSON serialization for embedded 2D direct lattices.

The wire format stores only the defining ambient :math:`3\times2` basis and its
unit. Induced metrics, dual bases, frames, and projectors are deliberately not
serialized; they are reconstructed through
:class:`~projectkoios.physkit.periodic.lattice.embedding.EmbeddedDirectLattice2DAdapter`
from the validated defining basis.
"""

from __future__ import annotations

import json
from typing import ClassVar

import numpy as np

from projectkoios.physkit.periodic.lattice.embedding import EmbeddedDirectLattice2D
from projectkoios.physkit.units import MatrixQuantity, PhysicalUnit, Unitless

_SCHEMA_VERSION = 1
_REPRESENTATION = "embedded-direct-lattice-2d"
_VECTOR_AXIS = "columns"
_MAX_CONTENT_CHARACTERS = 16_384
_MAX_UNIT_EXPRESSION_CHARACTERS = 256


class EmbeddedDirectLattice2DSerializationError(ValueError):
    """Report a strict embedded-lattice JSON wire-contract violation."""


class EmbeddedDirectLattice2DJsonCodec:
    r"""Encode and decode the strict schema-version-1 embedded-plane record.

    The deterministic JSON object contains exactly ``schema_version``,
    ``representation``, and ``ambient_basis``. The basis object contains exactly
    ``a1``, ``a2``, ``vector_axis``, and ``unit``. Primitive vectors are
    three-number arrays and are explicitly declared as columns.

    Physical units use
    ``{"kind":"physical","expression":"..."}``; nondimensional bases use
    ``{"kind":"unitless"}``. Unknown, missing, or duplicate object members are
    rejected at every level. Numbers must be standard finite JSON numbers;
    booleans and numeric strings are not accepted as numbers.

    Attributes
    ----------
    SCHEMA_VERSION
        Only supported integer schema version, equal to ``1``.
    REPRESENTATION
        Exact representation discriminator ``"embedded-direct-lattice-2d"``.
    MAX_CONTENT_CHARACTERS
        Maximum accepted JSON text length. The fixed six-number schema needs no
        unbounded wire allocation.
    """

    __slots__ = ()

    SCHEMA_VERSION: ClassVar[int] = _SCHEMA_VERSION
    REPRESENTATION: ClassVar[str] = _REPRESENTATION
    MAX_CONTENT_CHARACTERS: ClassVar[int] = _MAX_CONTENT_CHARACTERS

    _TOP_LEVEL_FIELDS: ClassVar[frozenset[str]] = frozenset(
        {"schema_version", "representation", "ambient_basis"}
    )
    _BASIS_FIELDS: ClassVar[frozenset[str]] = frozenset(
        {"a1", "a2", "vector_axis", "unit"}
    )
    _PHYSICAL_UNIT_FIELDS: ClassVar[frozenset[str]] = frozenset({"kind", "expression"})
    _UNITLESS_FIELDS: ClassVar[frozenset[str]] = frozenset({"kind"})

    def dumps(self, embedded_lattice: EmbeddedDirectLattice2D) -> str:
        """Return deterministic schema-version-1 JSON text.

        Parameters
        ----------
        embedded_lattice
            Exact nominal embedded-lattice record to serialize.

        Returns
        -------
        str
            Compact sorted JSON ending in one newline. Only the defining basis
            and unit are emitted.

        Raises
        ------
        EmbeddedDirectLattice2DSerializationError
            If the input is not exactly `EmbeddedDirectLattice2D` or finite JSON
            cannot be emitted.
        """
        if type(embedded_lattice) is not EmbeddedDirectLattice2D:
            raise EmbeddedDirectLattice2DSerializationError(
                "embedded_lattice must be exactly EmbeddedDirectLattice2D"
            )

        basis = embedded_lattice.ambient_basis
        unit_payload: dict[str, object]
        if isinstance(basis.unit, Unitless):
            unit_payload = {"kind": "unitless"}
        else:
            # EmbeddedDirectLattice2D already restricts this branch to a
            # physical length unit; emit its exact Pint expression.
            assert isinstance(basis.unit, PhysicalUnit)
            unit_payload = {
                "expression": basis.unit.expression,
                "kind": "physical",
            }
        payload = {
            "ambient_basis": {
                "a1": basis.magnitude[:, 0].tolist(),
                "a2": basis.magnitude[:, 1].tolist(),
                "unit": unit_payload,
                "vector_axis": _VECTOR_AXIS,
            },
            "representation": _REPRESENTATION,
            "schema_version": _SCHEMA_VERSION,
        }
        try:
            content = (
                json.dumps(
                    payload,
                    sort_keys=True,
                    separators=(",", ":"),
                    allow_nan=False,
                )
                + "\n"
            )
        except (TypeError, ValueError) as error:
            raise EmbeddedDirectLattice2DSerializationError(
                "embedded lattice cannot be represented as finite JSON"
            ) from error
        if len(content) > _MAX_CONTENT_CHARACTERS:
            raise EmbeddedDirectLattice2DSerializationError(
                "embedded-lattice JSON exceeds the resource bound"
            )
        return content

    def loads(self, content: str) -> EmbeddedDirectLattice2D:
        """Reconstruct an exact embedded-lattice record from strict JSON.

        Parameters
        ----------
        content
            Schema-version-1 JSON text. Its character count must not exceed
            `MAX_CONTENT_CHARACTERS`.

        Returns
        -------
        EmbeddedDirectLattice2D
            Newly reconstructed exact nominal record. Derived geometry is not
            accepted from the wire and can be recomputed by the adapter.

        Raises
        ------
        EmbeddedDirectLattice2DSerializationError
            If text is malformed or oversized; contains duplicate, unknown, or
            missing members; uses unsupported schema metadata; contains invalid
            numbers, vectors, units, or basis geometry; or otherwise violates
            the domain record invariants.
        """
        if type(content) is not str:
            raise EmbeddedDirectLattice2DSerializationError("content must be a string")
        if len(content) > _MAX_CONTENT_CHARACTERS:
            raise EmbeddedDirectLattice2DSerializationError(
                "embedded-lattice JSON exceeds the resource bound"
            )
        try:
            payload = json.loads(
                content,
                object_pairs_hook=self._unique_object,
                parse_constant=self._reject_constant,
                parse_int=self._parse_integer,
            )
        except EmbeddedDirectLattice2DSerializationError:
            raise
        except (json.JSONDecodeError, RecursionError, ValueError) as error:
            raise EmbeddedDirectLattice2DSerializationError(
                "embedded-lattice content is not strict valid JSON"
            ) from error
        if type(payload) is not dict:
            raise EmbeddedDirectLattice2DSerializationError(
                "embedded-lattice document must be an object"
            )

        self._require_fields(payload, self._TOP_LEVEL_FIELDS, "document")
        if self._require_integer(payload["schema_version"], "schema_version") != 1:
            raise EmbeddedDirectLattice2DSerializationError(
                "unsupported embedded-lattice JSON schema"
            )
        if (
            self._require_string(payload["representation"], "representation")
            != _REPRESENTATION
        ):
            raise EmbeddedDirectLattice2DSerializationError(
                "unsupported embedded-lattice representation"
            )

        basis = self._require_object(payload["ambient_basis"], "ambient_basis")
        self._require_fields(basis, self._BASIS_FIELDS, "ambient_basis")
        if (
            self._require_string(
                basis["vector_axis"],
                "ambient_basis.vector_axis",
            )
            != _VECTOR_AXIS
        ):
            raise EmbeddedDirectLattice2DSerializationError(
                "ambient_basis.vector_axis must declare vectors as columns"
            )
        unit = self._decode_unit(basis["unit"])
        a1 = self._decode_vector(basis["a1"], "ambient_basis.a1")
        a2 = self._decode_vector(basis["a2"], "ambient_basis.a2")

        # Reconstruct only the defining basis. Domain constructors reapply
        # finiteness, unit, shape, and scale-invariant rank invariants.
        try:
            return EmbeddedDirectLattice2D(
                ambient_basis=MatrixQuantity(
                    np.column_stack((a1, a2)),
                    unit,
                )
            )
        except (TypeError, ValueError) as error:
            raise EmbeddedDirectLattice2DSerializationError(
                "decoded values do not form a valid embedded lattice"
            ) from error

    def _decode_unit(self, value: object) -> PhysicalUnit | Unitless:
        """Decode one exact discriminated basis-unit object."""
        unit = self._require_object(value, "ambient_basis.unit")
        if "kind" not in unit:
            raise EmbeddedDirectLattice2DSerializationError(
                "ambient_basis.unit fields are invalid: missing ['kind']"
            )
        kind = self._require_string(unit["kind"], "ambient_basis.unit.kind")
        if kind == "unitless":
            self._require_fields(
                unit,
                self._UNITLESS_FIELDS,
                "ambient_basis.unit",
            )
            return Unitless()
        if kind == "physical":
            self._require_fields(
                unit,
                self._PHYSICAL_UNIT_FIELDS,
                "ambient_basis.unit",
            )
            expression = self._require_string(
                unit["expression"],
                "ambient_basis.unit.expression",
            )
            if len(expression) > _MAX_UNIT_EXPRESSION_CHARACTERS:
                raise EmbeddedDirectLattice2DSerializationError(
                    "ambient_basis.unit.expression exceeds the resource bound"
                )
            try:
                return PhysicalUnit(expression)
            except (TypeError, ValueError) as error:
                raise EmbeddedDirectLattice2DSerializationError(
                    "ambient_basis.unit.expression must be a physical unit"
                ) from error
        raise EmbeddedDirectLattice2DSerializationError(
            "ambient_basis.unit.kind must be 'physical' or 'unitless'"
        )

    def _decode_vector(self, value: object, path: str) -> np.ndarray:
        """Decode one exact three-number ambient primitive vector."""
        if type(value) is not list or len(value) != 3:
            raise EmbeddedDirectLattice2DSerializationError(
                f"{path} must contain exactly three finite JSON numbers"
            )
        return np.asarray(
            [
                self._require_real(component, f"{path}[{index}]")
                for index, component in enumerate(value)
            ],
            dtype=np.float64,
        )

    @staticmethod
    def _require_fields(
        value: dict[str, object],
        required: frozenset[str],
        path: str,
    ) -> None:
        """Require one JSON object to have exactly the declared member names."""
        fields = frozenset(value)
        missing = sorted(required - fields)
        unknown = sorted(fields - required)
        if missing or unknown:
            details: list[str] = []
            if missing:
                details.append(f"missing {missing}")
            if unknown:
                details.append(f"unknown {unknown}")
            raise EmbeddedDirectLattice2DSerializationError(
                f"{path} fields are invalid: {', '.join(details)}"
            )

    @staticmethod
    def _require_object(value: object, path: str) -> dict[str, object]:
        """Return one exact JSON object."""
        if type(value) is not dict:
            raise EmbeddedDirectLattice2DSerializationError(f"{path} must be an object")
        return value

    @staticmethod
    def _require_integer(value: object, path: str) -> int:
        """Return an exact JSON integer while rejecting booleans and strings."""
        if type(value) is not int:
            raise EmbeddedDirectLattice2DSerializationError(
                f"{path} must be a JSON integer"
            )
        return value

    @staticmethod
    def _require_real(value: object, path: str) -> float:
        """Return one finite JSON real while rejecting booleans and strings."""
        if type(value) not in (int, float):
            raise EmbeddedDirectLattice2DSerializationError(
                f"{path} must be a JSON number"
            )
        try:
            number = float(value)
        except OverflowError as error:
            raise EmbeddedDirectLattice2DSerializationError(
                f"{path} must be finite"
            ) from error
        if not np.isfinite(number):
            raise EmbeddedDirectLattice2DSerializationError(f"{path} must be finite")
        return number

    @staticmethod
    def _require_string(value: object, path: str) -> str:
        """Return one nonempty stripped JSON string."""
        if type(value) is not str or not value or value != value.strip():
            raise EmbeddedDirectLattice2DSerializationError(
                f"{path} must be a nonempty stripped string"
            )
        return value

    @staticmethod
    def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
        """Reject duplicate JSON member names before dictionary collapse."""
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise EmbeddedDirectLattice2DSerializationError(
                    f"duplicate JSON object key: {key!r}"
                )
            result[key] = value
        return result

    @staticmethod
    def _reject_constant(constant: str) -> None:
        """Reject nonstandard ``NaN`` and infinity constants."""
        raise EmbeddedDirectLattice2DSerializationError(
            f"nonstandard JSON constant is not allowed: {constant}"
        )

    @staticmethod
    def _parse_integer(token: str) -> int:
        """Parse one bounded-content JSON integer token."""
        try:
            return int(token)
        except ValueError as error:
            raise EmbeddedDirectLattice2DSerializationError(
                "JSON integer token is invalid or exceeds the resource bound"
            ) from error
