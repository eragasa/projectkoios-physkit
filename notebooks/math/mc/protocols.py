from typing import Protocol


class State(Protocol):
    """One configuration of a system."""

    ...


class StateSpace(Protocol):
    """Collection of valid states for a model."""

    def contains(self, state: State) -> bool:
        """Return True if state belongs to this state space."""
        ...
