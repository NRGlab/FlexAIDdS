"""flexaidds: Python bindings and read-only analysis helpers for FlexAID∆S."""

from .thermodynamics import StatMechEngine, Thermodynamics
from .models import BindingModeResult, DockingResult, PoseResult
from .results import load_results

__all__ = [
    "StatMechEngine",
    "Thermodynamics",
    "PoseResult",
    "BindingModeResult",
    "DockingResult",
    "load_results",
]
__version__ = "0.1.0"
