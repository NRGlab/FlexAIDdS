"""High-level docking interface for FlexAID∆S.

Provides Pythonic API for molecular docking workflows.
"""

import numpy as np
from pathlib import Path
from typing import List, Optional, Dict, Any
from dataclasses import dataclass

from .thermodynamics import Thermodynamics, StatMechEngine

try:
    from . import _core
except ImportError:
    _core = None


@dataclass
class Pose:
    """Single docked pose within a binding mode.
    
    Attributes:
        index: Pose index in GA population
        energy: Binding energy (CF score) in kcal/mol
        rmsd: RMSD to reference structure (if available)
        coordinates: Atomic coordinates (Nx3 array)
        boltzmann_weight: Statistical weight in ensemble
    """
    index: int
    energy: float
    rmsd: Optional[float] = None
    coordinates: Optional[np.ndarray] = None
    boltzmann_weight: float = 0.0
    
    def to_dict(self) -> dict:
        return {
            'index': self.index,
            'energy_kcal_mol': self.energy,
            'rmsd_angstrom': self.rmsd,
            'boltzmann_weight': self.boltzmann_weight,
        }


class BindingMode:
    """Binding mode: cluster of docked poses with thermodynamic scoring.
    
    A binding mode represents a distinct local minimum on the binding energy
    landscape, characterized by an ensemble of similar poses.
    
    Example:
        >>> mode = results.binding_modes[0]  # top-ranked mode
        >>> thermo = mode.get_thermodynamics()
        >>> print(f"ΔG = {thermo.free_energy:.2f} kcal/mol")
        >>> print(f"ΔH = {thermo.mean_energy:.2f}, TΔS = {thermo.entropy_term:.2f}")
    """
    
    def __init__(self, cpp_binding_mode=None):
        """Initialize from C++ BindingMode object (internal use)."""
        self._cpp_mode = cpp_binding_mode
        self._poses: List[Pose] = []
    
    def get_thermodynamics(self) -> Thermodynamics:
        """Get full thermodynamic properties of this binding mode.
        
        Returns:
            Thermodynamics object with F, S, H, Cv, etc.
        """
        if self._cpp_mode is None:
            raise RuntimeError("C++ binding mode not initialized")
        thermo_cpp = self._cpp_mode.get_thermodynamics()
        return Thermodynamics(
            temperature=thermo_cpp.temperature,
            log_Z=thermo_cpp.log_Z,
            free_energy=thermo_cpp.free_energy,
            mean_energy=thermo_cpp.mean_energy,
            mean_energy_sq=thermo_cpp.mean_energy_sq,
            heat_capacity=thermo_cpp.heat_capacity,
            entropy=thermo_cpp.entropy,
            std_energy=thermo_cpp.std_energy,
        )
    
    @property
    def free_energy(self) -> float:
        """Helmholtz free energy F = H - TS (kcal/mol)."""
        if self._cpp_mode:
            return self._cpp_mode.get_free_energy()
        return float('inf')
    
    @property
    def enthalpy(self) -> float:
        """Boltzmann-weighted average energy ⟨E⟩ (kcal/mol)."""
        if self._cpp_mode:
            return self._cpp_mode.compute_enthalpy()
        return float('inf')
    
    @property
    def entropy(self) -> float:
        """Configurational entropy S (kcal mol⁻¹ K⁻¹)."""
        if self._cpp_mode:
            return self._cpp_mode.compute_entropy()
        return 0.0
    
    @property
    def n_poses(self) -> int:
        """Number of poses in this binding mode."""
        if self._cpp_mode:
            return self._cpp_mode.get_BindingMode_size()
        return len(self._poses)
    
    def __len__(self) -> int:
        return self.n_poses
    
    def __repr__(self) -> str:
        return (f"<BindingMode n_poses={self.n_poses} "
                f"F={self.free_energy:.2f} H={self.enthalpy:.2f} "
                f"S={self.entropy:.5f}>")


class BindingPopulation:
    """Collection of binding modes from a docking run.
    
    Provides ensemble-level analysis and ranking of binding modes.
    """
    
    def __init__(self):
        self._modes: List[BindingMode] = []
        self._temperature: float = 300.0
    
    def add_mode(self, mode: BindingMode) -> None:
        """Add a binding mode to the population."""
        self._modes.append(mode)
    
    def rank_by_free_energy(self) -> List[BindingMode]:
        """Return binding modes sorted by free energy (best first)."""
        return sorted(self._modes, key=lambda m: m.free_energy)
    
    def compute_global_thermodynamics(self) -> Thermodynamics:
        """Compute thermodynamics over all binding modes.
        
        Returns:
            Global ensemble thermodynamics
        """
        engine = StatMechEngine(self._temperature)
        for mode in self._modes:
            # Aggregate all pose energies from all modes
            for _ in range(mode.n_poses):
                engine.add_sample(mode.enthalpy)  # Simplified: use mode average
        return engine.compute()
    
    @property
    def n_modes(self) -> int:
        """Number of binding modes."""
        return len(self._modes)
    
    def __len__(self) -> int:
        return self.n_modes
    
    def __getitem__(self, index: int) -> BindingMode:
        return self._modes[index]
    
    def __iter__(self):
        return iter(self._modes)
    
    def __repr__(self) -> str:
        return f"<BindingPopulation n_modes={self.n_modes} T={self._temperature}K>"


class Docking:
    """High-level interface for FlexAID∆S molecular docking.
    
    Example:
        >>> docking = Docking("config.inp")
        >>> results = docking.run()
        >>> top_mode = results.binding_modes[0]
        >>> print(f"Best ΔG: {top_mode.free_energy:.2f} kcal/mol")
    """
    
    def __init__(self, config_file: str):
        """Initialize docking from configuration file.
        
        Args:
            config_file: Path to FlexAID .inp config file
        """
        self.config_file = Path(config_file)
        if not self.config_file.exists():
            raise FileNotFoundError(f"Config file not found: {config_file}")
        
        self._config: Dict[str, Any] = {}
        self._parse_config()
    
    def _parse_config(self) -> None:
        """Parse FlexAID config file (stub - needs implementation)."""
        # TODO: Implement config parser
        # This will read INPLIG, PDBNAM, RNGOPT, OPTIMZ, etc.
        pass
    
    def run(self, **kwargs) -> BindingPopulation:
        """Execute docking simulation.
        
        Returns:
            BindingPopulation with ranked binding modes
        
        Note:
            Full implementation requires integration with C++ FlexAID GA engine.
            This is a stub for Phase 2 development.
        """
        raise NotImplementedError(
            "Full docking pipeline integration in progress. "
            "For now, use C++ FlexAID binary directly and parse output with "
            "BindingMode/BindingPopulation wrappers."
        )
    
    def __repr__(self) -> str:
        return f"<Docking config={self.config_file.name}>"
