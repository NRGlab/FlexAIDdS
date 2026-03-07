"""PyMOL visualization functions for FlexAID∆S binding modes.

These functions can be called from PyMOL command line:
    PyMOL> flexaids_load /path/to/output
    PyMOL> flexaids_show_ensemble mode1
    PyMOL> flexaids_color_boltzmann mode1
    PyMOL> flexaids_thermo mode1
"""

import os
import re
import math
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple

try:
    from pymol import cmd, stored
    import pymol
except ImportError:
    raise ImportError("PyMOL not available")


# ─── physical constant ────────────────────────────────────────────────────────
_kB_kcal: float = 0.001987206  # kcal mol⁻¹ K⁻¹


# ─── data model ──────────────────────────────────────────────────────────────
@dataclass
class _ModeRecord:
    """Thermodynamic and structural data for one FlexAID binding mode."""
    mode_id: int
    pdb_objects: List[str] = field(default_factory=list)  # PyMOL object names
    cf_values: List[float] = field(default_factory=list)  # CF score per pose
    best_cf: float = 0.0
    total_cf: float = 0.0
    frequency: int = 0
    # Thermodynamics (populated after parsing)
    free_energy: Optional[float] = None   # F = -kT ln Z  (kcal/mol)
    enthalpy: Optional[float] = None      # <E>            (kcal/mol)
    entropy: Optional[float] = None       # S              (kcal mol⁻¹ K⁻¹)
    heat_capacity: Optional[float] = None # Cv             (kcal mol⁻¹ K⁻²)
    boltzmann_weights: List[float] = field(default_factory=list)


# Global state: keyed by mode_name (e.g. "mode1")
_loaded_modes: Dict[str, _ModeRecord] = {}
_output_dir: Optional[Path] = None
_temperature_K: float = 300.0


# ─── private helpers ─────────────────────────────────────────────────────────

def _parse_pdb_cf_values(pdb_path: Path) -> Tuple[Optional[int], List[float], float, float, int]:
    """Extract CF values and mode metadata from a FlexAID result PDB.

    Returns:
        (mode_id, cf_values, best_cf, total_cf, frequency)
        mode_id may be None if the REMARK is absent.
    """
    mode_id: Optional[int] = None
    cf_values: List[float] = []
    best_cf: float = 0.0
    total_cf: float = 0.0
    frequency: int = 0

    _re_mode = re.compile(
        r"REMARK\s+Binding Mode:(\d+)"
        r"\s+Best CF in Binding Mode:([\d.+-]+)"
        r".*?Binding Mode Total CF:([\d.+-]+)"
        r"\s+Binding Mode Frequency:(\d+)",
    )
    _re_cf = re.compile(r"REMARK\s+CF=([\d.eE+\-]+)")

    try:
        with open(pdb_path) as fh:
            for line in fh:
                m = _re_mode.search(line)
                if m:
                    mode_id = int(m.group(1))
                    best_cf = float(m.group(2))
                    total_cf = float(m.group(3))
                    frequency = int(m.group(4))
                    continue
                m = _re_cf.search(line)
                if m:
                    cf_values.append(float(m.group(1)))
    except OSError:
        pass

    # Static (non-dynamic) result: single representative pose, one CF in REMARK
    if not cf_values and best_cf != 0.0:
        cf_values = [best_cf]

    return mode_id, cf_values, best_cf, total_cf, frequency


def _parse_cad_file(cad_path: Path) -> Dict[int, Dict[str, Any]]:
    """Parse a FlexAID .cad cluster-assignment file.

    Returns:
        {cluster_id: {'best_cf': float, 'total_cf': float, 'frequency': int}}
    """
    result: Dict[int, Dict[str, Any]] = {}
    _re = re.compile(
        r"Cluster\s+(\d+):\s+Center:\d+\s+\(CF:([\d.eE+\-]+)\)"
        r"\s+Best:\d+\s+\(CF:([\d.eE+\-]+)\)"
        r"\s+TCF=([\d.eE+\-]+)"
        r"\s+Frequency=(\d+)",
    )
    try:
        with open(cad_path) as fh:
            for line in fh:
                m = _re.search(line)
                if m:
                    cid = int(m.group(1))
                    result[cid] = {
                        "best_cf": float(m.group(3)),
                        "total_cf": float(m.group(4)),
                        "frequency": int(m.group(5)),
                    }
    except OSError:
        pass
    return result


def _compute_thermodynamics(cf_values: List[float], temperature_K: float) -> Dict[str, float]:
    """Pure-Python canonical-ensemble thermodynamics from CF-score samples.

    Uses numerically stable log-sum-exp to avoid overflow/underflow.

    Returns dict with keys: free_energy, enthalpy, entropy, heat_capacity,
    boltzmann_weights (list).
    """
    if not cf_values:
        return {}

    beta = 1.0 / (_kB_kcal * temperature_K)
    energies = cf_values

    # Numerically stable partition function via log-sum-exp
    neg_beta_E = [-beta * e for e in energies]
    max_val = max(neg_beta_E)
    shifted = [math.exp(v - max_val) for v in neg_beta_E]
    sum_shifted = sum(shifted)
    log_Z = max_val + math.log(sum_shifted)

    # Boltzmann weights
    weights = [s / sum_shifted for s in shifted]

    # Mean energy <E>
    mean_E = sum(w * e for w, e in zip(weights, energies))

    # Mean energy squared <E²>
    mean_E2 = sum(w * e * e for w, e in zip(weights, energies))

    # Thermodynamic quantities
    free_energy = -_kB_kcal * temperature_K * log_Z
    entropy = (mean_E - free_energy) / temperature_K
    heat_capacity = (mean_E2 - mean_E ** 2) / (_kB_kcal * temperature_K ** 2)

    return {
        "free_energy": free_energy,
        "enthalpy": mean_E,
        "entropy": entropy,
        "heat_capacity": heat_capacity,
        "boltzmann_weights": weights,
    }


# ─── public API ──────────────────────────────────────────────────────────────

def load_binding_modes(output_dir: str, temperature: float = 300.0) -> None:
    """Load FlexAID∆S docking results from output directory.

    Parses all ``result_*.pdb`` files and the accompanying ``.cad`` cluster
    file, computes canonical-ensemble thermodynamics for each binding mode,
    and loads structures into PyMOL.

    Args:
        output_dir: Path to FlexAID output directory
        temperature: Simulation temperature in Kelvin (default 300 K)

    Example:
        PyMOL> flexaids_load /data/docking_results
        PyMOL> flexaids_load /data/docking_results, 310
    """
    global _loaded_modes, _output_dir, _temperature_K

    output_path = Path(output_dir)
    if not output_path.exists():
        print(f"ERROR: Directory not found: {output_dir}")
        return

    _output_dir = output_path
    _temperature_K = float(temperature)
    _loaded_modes.clear()

    # ── locate result PDB files (exclude receptor/ligand .inp.pdb files) ──
    pdb_files = sorted(
        p for p in output_path.glob("result_*.pdb")
        if not p.name.endswith(".inp.pdb")
    )
    if not pdb_files:
        print(f"ERROR: No result_*.pdb files found in {output_dir}")
        return

    # ── parse optional .cad cluster metadata ──
    cad_files = list(output_path.glob("*.cad"))
    cad_meta: Dict[int, Dict[str, Any]] = {}
    if cad_files:
        cad_meta = _parse_cad_file(cad_files[0])

    # ── parse each PDB, group into _ModeRecord objects ──
    for pdb_file in pdb_files:
        mode_id, cf_values, best_cf, total_cf, frequency = _parse_pdb_cf_values(pdb_file)

        # Fall back to sequential numbering if REMARK absent
        if mode_id is None:
            stem_parts = pdb_file.stem.rsplit("_", 1)
            try:
                mode_id = int(stem_parts[-1])
            except ValueError:
                mode_id = len(_loaded_modes) + 1

        mode_name = f"mode{mode_id}"
        if mode_name not in _loaded_modes:
            cad = cad_meta.get(mode_id, {})
            _loaded_modes[mode_name] = _ModeRecord(
                mode_id=mode_id,
                best_cf=cad.get("best_cf", best_cf),
                total_cf=cad.get("total_cf", total_cf),
                frequency=cad.get("frequency", frequency),
            )

        rec = _loaded_modes[mode_name]
        rec.cf_values.extend(cf_values)

        # Load into PyMOL
        obj_name = f"flexaids_{mode_name}_{pdb_file.stem}"
        cmd.load(str(pdb_file), obj_name)
        cmd.disable(obj_name)
        rec.pdb_objects.append(obj_name)

    # ── compute thermodynamics per mode ──
    for mode_name, rec in _loaded_modes.items():
        if rec.cf_values:
            td = _compute_thermodynamics(rec.cf_values, _temperature_K)
            rec.free_energy = td.get("free_energy")
            rec.enthalpy = td.get("enthalpy")
            rec.entropy = td.get("entropy")
            rec.heat_capacity = td.get("heat_capacity")
            rec.boltzmann_weights = td.get("boltzmann_weights", [])
            rec.frequency = len(rec.cf_values)

    n_modes = len(_loaded_modes)
    n_poses = sum(len(r.pdb_objects) for r in _loaded_modes.values())
    print(f"Loaded {n_modes} binding modes ({n_poses} PDB objects) from {output_dir}")
    print("Use 'flexaids_show_ensemble modeN' to visualize a binding mode.")


def show_pose_ensemble(mode_name: str, show_all: bool = True) -> None:
    """Display all poses belonging to a binding mode.

    Args:
        mode_name: Binding mode identifier (e.g. 'mode1')
        show_all: If True, show all poses; if False, show representative only
            (the pose with the highest Boltzmann weight).

    Example:
        PyMOL> flexaids_show_ensemble mode1
        PyMOL> flexaids_show_ensemble mode1, show_all=0
    """
    if not _loaded_modes:
        print("ERROR: No modes loaded. Use 'flexaids_load' first.")
        return

    rec = _loaded_modes.get(mode_name)
    if rec is None:
        available = ", ".join(sorted(_loaded_modes))
        print(f"ERROR: Mode '{mode_name}' not found. Available: {available}")
        return

    if not rec.pdb_objects:
        print(f"ERROR: No PDB objects for {mode_name}.")
        return

    if show_all:
        for obj in rec.pdb_objects:
            cmd.enable(obj)
            cmd.show("cartoon", obj)
            cmd.show("sticks", f"{obj} and organic")
    else:
        # Representative = object corresponding to highest Boltzmann weight.
        # Weights are indexed parallel to pdb_objects (one weight per file).
        if rec.boltzmann_weights and len(rec.boltzmann_weights) == len(rec.pdb_objects):
            rep_idx = rec.boltzmann_weights.index(max(rec.boltzmann_weights))
        else:
            rep_idx = 0  # fall back to first (lowest CF)

        for i, obj in enumerate(rec.pdb_objects):
            if i == rep_idx:
                cmd.enable(obj)
                cmd.show("cartoon", obj)
                cmd.show("sticks", f"{obj} and organic")
            else:
                cmd.disable(obj)

    zoom_sel = " ".join(rec.pdb_objects)
    cmd.zoom(zoom_sel)
    label = "all poses" if show_all else "representative pose"
    print(f"Showing {label} for {mode_name} ({len(rec.pdb_objects)} PDB objects).")


def color_by_boltzmann_weight(mode_name: str) -> None:
    """Color poses by Boltzmann weight (blue = low probability, red = high).

    Args:
        mode_name: Binding mode identifier

    Example:
        PyMOL> flexaids_color_boltzmann mode1
    """
    if not _loaded_modes:
        print("ERROR: No modes loaded. Use 'flexaids_load' first.")
        return

    rec = _loaded_modes.get(mode_name)
    if rec is None:
        available = ", ".join(sorted(_loaded_modes))
        print(f"ERROR: Mode '{mode_name}' not found. Available: {available}")
        return

    pose_objects = rec.pdb_objects
    if not pose_objects:
        print(f"ERROR: No poses for {mode_name}.")
        return

    # Use computed Boltzmann weights if available, otherwise fall back to
    # uniform weights (all poses equally probable).
    weights = rec.boltzmann_weights
    if not weights or len(weights) != len(pose_objects):
        n = len(pose_objects)
        weights = [1.0 / n] * n

    w_min = min(weights)
    w_max = max(weights)
    w_range = w_max - w_min if w_max > w_min else 1.0

    for i, (obj, w) in enumerate(zip(pose_objects, weights)):
        t = (w - w_min) / w_range  # 0 = lowest weight (blue), 1 = highest (red)
        color_name = f"flexaids_bw_{mode_name}_{i}"
        cmd.set_color(color_name, [t, 0.0, 1.0 - t])
        cmd.color(color_name, obj)
        cmd.enable(obj)

    print(f"Colored {len(pose_objects)} poses for {mode_name} by Boltzmann weight "
          f"(blue=low, red=high).")


def show_thermodynamics(mode_name: str) -> None:
    """Print thermodynamic properties of a binding mode to PyMOL console.

    Args:
        mode_name: Binding mode identifier

    Example:
        PyMOL> flexaids_thermo mode1
    """
    if not _loaded_modes:
        print("ERROR: No modes loaded. Use 'flexaids_load' first.")
        return

    rec = _loaded_modes.get(mode_name)
    if rec is None:
        available = ", ".join(sorted(_loaded_modes))
        print(f"ERROR: Mode '{mode_name}' not found. Available: {available}")
        return

    T = _temperature_K
    entropy_term = (rec.entropy * T) if rec.entropy is not None else None

    print(f"\nThermodynamics for {mode_name} (T = {T:.1f} K):")
    print(f"  ΔG (Free Energy):     {rec.free_energy:10.4f} kcal/mol"
          if rec.free_energy is not None else "  ΔG (Free Energy):     N/A")
    print(f"  ΔH (Enthalpy):        {rec.enthalpy:10.4f} kcal/mol"
          if rec.enthalpy is not None else "  ΔH (Enthalpy):        N/A")
    print(f"  S (Entropy):          {rec.entropy:10.6f} kcal/(mol·K)"
          if rec.entropy is not None else "  S (Entropy):          N/A")
    print(f"  TΔS (Entropy term):   {entropy_term:10.4f} kcal/mol"
          if entropy_term is not None else "  TΔS (Entropy term):   N/A")
    print(f"  Heat Capacity (Cv):   {rec.heat_capacity:10.4f} kcal/(mol·K²)"
          if rec.heat_capacity is not None else "  Heat Capacity (Cv):   N/A")
    print(f"  Best CF score:        {rec.best_cf:10.5f}")
    print(f"  # Poses:              {rec.frequency:10d}")
    print()


def export_to_nrgsuite(output_dir: str, nrgsuite_file: str) -> None:
    """Export binding modes to NRGSuite-compatible format.

    Writes a tab-separated text file with one row per binding mode containing
    the mode ID, representative CF score, free energy, enthalpy, entropy, and
    pose count.  This format is compatible with the NRGSuite FreeNRG pipeline.

    Args:
        output_dir: FlexAID output directory (used to resolve PDB paths)
        nrgsuite_file: Output file path for NRGSuite

    Example:
        PyMOL> flexaids_export /data/docking_results, /data/nrgsuite_input.txt
    """
    if not _loaded_modes:
        # Attempt lazy load if the caller passes a directory we haven't loaded.
        load_binding_modes(output_dir)
        if not _loaded_modes:
            print(f"ERROR: Could not load any modes from {output_dir}")
            return

    out_path = Path(nrgsuite_file)
    try:
        with open(out_path, "w") as fh:
            fh.write(
                "# FlexAID∆S → NRGSuite export\n"
                "# mode_id\tbest_cf\tfree_energy_kcal_mol\t"
                "enthalpy_kcal_mol\tentropy_kcal_mol_K\tn_poses\n"
            )
            for mode_name in sorted(_loaded_modes, key=lambda n: _loaded_modes[n].mode_id):
                rec = _loaded_modes[mode_name]
                fh.write(
                    f"{rec.mode_id}\t"
                    f"{rec.best_cf:.5f}\t"
                    f"{rec.free_energy:.4f}\t"
                    f"{rec.enthalpy:.4f}\t"
                    f"{rec.entropy:.6f}\t"
                    f"{rec.frequency}\n"
                    if None not in (rec.free_energy, rec.enthalpy, rec.entropy)
                    else f"{rec.mode_id}\t{rec.best_cf:.5f}\tN/A\tN/A\tN/A\t{rec.frequency}\n"
                )
    except OSError as exc:
        print(f"ERROR: Could not write NRGSuite file: {exc}")
        return

    print(f"Exported {len(_loaded_modes)} binding modes to {out_path}")


# ─── Auto-register PyMOL commands when module is imported ────────────────────
cmd.extend("flexaids_load", load_binding_modes)
cmd.extend("flexaids_show_ensemble", show_pose_ensemble)
cmd.extend("flexaids_color_boltzmann", color_by_boltzmann_weight)
cmd.extend("flexaids_thermo", show_thermodynamics)
cmd.extend("flexaids_export", export_to_nrgsuite)
