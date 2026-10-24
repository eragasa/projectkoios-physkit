"""Explicit failures for unavailable represented spectral decompositions."""


class UnsupportedSpectralDecompositionError(NotImplementedError):
    """Report a valid matrix family without an implemented decomposition."""
