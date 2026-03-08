"""flexaidds: Python bindings and read-only analysis helpers for FlexAID∆S."""

try:
    from ._core import StatMechEngine, Thermodynamics
except ImportError:
    StatMechEngine = None  # type: ignore[assignment,misc]
    Thermodynamics = None  # type: ignore[assignment,misc]

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
