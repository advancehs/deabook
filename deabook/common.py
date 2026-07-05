"""Public common helpers for object-oriented DEA implementations."""

from .core import (
    CNLSResultMixin,
    ModelDictionaryMixin,
    PyomoResultMixin,
    format_ref_indexed_results,
    format_var_indexed_results,
    indexed_component_to_frame,
    indexed_component_to_numpy,
    scalar_component_to_series,
)


def hello_world():
    """Prints ``"Hello World!"`` to the console."""
    print("Hello World!")


__all__ = [
    "CNLSResultMixin",
    "ModelDictionaryMixin",
    "PyomoResultMixin",
    "format_ref_indexed_results",
    "format_var_indexed_results",
    "indexed_component_to_frame",
    "indexed_component_to_numpy",
    "scalar_component_to_series",
    "hello_world",
]
