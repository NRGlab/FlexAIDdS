"""flexaidds: Python bindings and read-only analysis helpers for FlexAID∆S."""

from .__version__ import __version__
from .thermodynamics import StatMechEngine, Thermodynamics
from .models import BindingModeResult, DockingResult, PoseResult
from .results import load_results

__all__ = [
    "__version__",
    "StatMechEngine",
    "Thermodynamics",
    "PoseResult",
    "BindingModeResult",
    "DockingResult",
    "load_results",
]
