"""FlexAID∆S: Entropy-driven molecular docking with flexible side-chains.

A modernized C++20 implementation of FlexAID with:
- Statistical mechanics scoring (partition functions, free energies)
- Configurational and vibrational entropy corrections
- GPU acceleration (CUDA/Metal)
- Python bindings for high-level workflows

Example:
    >>> import flexaidds as fds
    >>> docking = fds.Docking("config.inp")
    >>> results = docking.run()
    >>> for mode in results.binding_modes:
    ...     print(f"ΔG = {mode.free_energy:.2f} kcal/mol")

Modules:
    thermodynamics: Statistical mechanics engine and free energy calculations
    docking: High-level docking interface
    visualization: PyMOL integration for pose rendering
    io: File I/O utilities (PDB, MOL2, config files)
"""

from .__version__ import __version__

try:
    from . import _core  # pybind11 compiled module
except ImportError as e:
    import warnings
    warnings.warn(
        f"Could not import C++ extension module: {e}. "
        "FlexAID∆S bindings not available. Build with 'pip install -e .'.",
        ImportWarning
    )
    _core = None

# High-level Python API
from .thermodynamics import (
    StatMechEngine,
    Thermodynamics,
    BoltzmannLUT,
)
from .docking import (
    Docking,
    BindingMode,
    BindingPopulation,
)

__all__ = [
    "__version__",
    # Thermodynamics
    "StatMechEngine",
    "Thermodynamics",
    "BoltzmannLUT",
    # Docking
    "Docking",
    "BindingMode",
    "BindingPopulation",
]
