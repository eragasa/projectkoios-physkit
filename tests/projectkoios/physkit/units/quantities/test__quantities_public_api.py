"""Public façade tests for canonical model-system unit identities."""

import projectkoios.physkit.units as units
import projectkoios.physkit.units.quantities as quantities


def test__units_facade__retains_quantity_owner_identities() -> None:
    """The package façade exports exact defining-module objects."""
    names = (
        "MODEL_SYSTEM_PINT_REGISTRY",
        "ModelSystemUnitConverter",
        "ComplexMatrixQuantity",
        "ComplexSparseMatrixQuantity",
        "ComplexVectorQuantity",
        "MatrixQuantity",
        "ModelSystemQuantity",
        "ModelSystemUnit",
        "PhysicalUnit",
        "PintUnitConverter",
        "ScalarQuantity",
        "SparseMatrixQuantity",
        "Unitless",
        "VectorQuantity",
    )

    for name in names:
        assert getattr(units, name) is getattr(quantities, name)
    assert quantities.__all__ == list(names)


def test__converter_instances__share_the_canonical_registry() -> None:
    """Ad hoc and shared converters use one registry identity."""
    converter = quantities.PintUnitConverter()

    assert converter.registry is quantities.MODEL_SYSTEM_PINT_REGISTRY
    assert (
        quantities.ModelSystemUnitConverter.registry
        is quantities.MODEL_SYSTEM_PINT_REGISTRY
    )
