from grid_1d import GridType1D


def test_members_exist():
    assert GridType1D.LEFT_CLOSED is not None
    assert GridType1D.RIGHT_CLOSED is not None
    assert GridType1D.OPEN is not None
    assert GridType1D.INTERIOR is not None
    assert GridType1D.MIDPOINT is not None
    assert GridType1D.CLOSED is not None


def test_is_enum_and_unique():
    # Enum members are unique identities
    assert GridType1D.LEFT_CLOSED is GridType1D.LEFT_CLOSED
    assert GridType1D.LEFT_CLOSED is not GridType1D.RIGHT_CLOSED


def test_iterates_in_definition_order():
    # Only if you care about order; Enum preserves definition order in iteration
    names = [m.name for m in GridType1D]
    assert names == [
        "LEFT_CLOSED",
        "RIGHT_CLOSED",
        "OPEN",
        "INTERIOR",
        "MIDPOINT",
        "CLOSED",
    ]
