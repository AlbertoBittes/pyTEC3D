"""APEC-M Python package."""

from importlib import import_module

from .apec_m import APEC_M
from .apec_m_model import APEC_M_Model
from .apec_algorithm import ApecAlgorithm

_PARAVIEW_EXPORTS = {
    "plot",
    "setup_paraview_console",
    "start_apec_for_paraview",
}

__all__ = [
    "APEC_M",
    "APEC_M_Model",
    "ApecAlgorithm",
    "plot",
    "setup_paraview_console",
    "start_apec_for_paraview",
]


def __getattr__(name: str):
    if name in _PARAVIEW_EXPORTS:
        module = import_module(".apec_m_paraview", __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals().keys()) | _PARAVIEW_EXPORTS)
