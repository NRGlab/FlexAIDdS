"""flexaidds: Python bindings and read-only analysis helpers for FlexAID∆S."""

try:
    from ._core import StatMechEngine as _CppStatMechEngine
    from ._core import Thermodynamics as _CppThermodynamics  # noqa: F401
except ImportError:
    _CppStatMechEngine = None  # type: ignore[assignment]
    _CppThermodynamics = None  # type: ignore[assignment]

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
