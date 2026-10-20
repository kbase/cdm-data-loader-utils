"""Shared error types for the converter packages."""


class ConversionError(ValueError):
    """Base class for all converter errors; direction-specific errors subclass this."""
