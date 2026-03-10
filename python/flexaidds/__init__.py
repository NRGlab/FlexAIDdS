"""flexaidds: Python bindings and read-only analysis helpers for FlexAID∆S."""

try:
    from ._core import (
        BoltzmannLUT,
        ENCoMEngine,
        NormalMode,
        Replica,
        State,
        StatMechEngine,
        Thermodynamics,
        TIPoint,
        VibrationalEntropy,
        WHAMBin,
        kB_kcal,
        kB_SI,
    )
    HAS_CORE_BINDINGS = True
except ImportError:
    StatMechEngine = None
    Thermodynamics = None
    State = None
    BoltzmannLUT = None
    Replica = None
    WHAMBin = None
    TIPoint = None
    ENCoMEngine = None
    NormalMode = None
    VibrationalEntropy = None
    kB_kcal = 0.0019872041
    kB_SI = 1.380649e-23
    HAS_CORE_BINDINGS = False

from .models import BindingModeResult, DockingResult, PoseResult
from .results import load_results

__all__ = [
    # C++ core: statistical mechanics
    "StatMechEngine",
    "Thermodynamics",
    "State",
    "BoltzmannLUT",
    # C++ core: parallel tempering & free energy methods
    "Replica",
    "WHAMBin",
    "TIPoint",
    # C++ core: ENCoM vibrational entropy
    "ENCoMEngine",
    "NormalMode",
    "VibrationalEntropy",
    # Physical constants
    "kB_kcal",
    "kB_SI",
    # Python models & I/O
    "PoseResult",
    "BindingModeResult",
    "DockingResult",
    "load_results",
    # Availability flag
    "HAS_CORE_BINDINGS",
]

__version__ = "0.1.0"
