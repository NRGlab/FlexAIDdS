from __future__ import annotations

import argparse
import json
from pathlib import Path

from .results import load_results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m flexaidds",
        description="Inspect FlexAID∆S docking result directories from Python.",
    )
    parser.add_argument("results_dir", type=Path, help="Directory containing docking result PDB files")
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable JSON instead of a human summary.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    result = load_results(args.results_dir)
    if args.json:
        payload = {
            "source_dir": str(result.source_dir),
            "temperature": result.temperature,
            "n_modes": result.n_modes,
            "metadata": result.metadata,
            "binding_modes": result.to_records(),
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    print(f"Results directory: {result.source_dir}")
    print(f"Binding modes: {result.n_modes}")
    if result.temperature is not None:
        print(f"Temperature: {result.temperature} K")
    print()

    if not result.binding_modes:
        return 0

    # Table header
    header = (
        f"{'mode':>5}  {'rank':>4}  {'poses':>5}  "
        f"{'best_cf':>10}  {'free_energy':>12}  "
        f"{'enthalpy':>10}  {'entropy':>12}"
    )
    print(header)
    print("-" * len(header))

    # Sort by rank for display
    top = result.top_mode()
    for mode in sorted(result.binding_modes, key=lambda m: m.rank):
        cf_str = f"{mode.best_cf:10.4f}" if mode.best_cf is not None else f"{'N/A':>10}"
        fe_str = f"{mode.free_energy:12.4f}" if mode.free_energy is not None else f"{'N/A':>12}"
        h_str = f"{mode.enthalpy:10.4f}" if mode.enthalpy is not None else f"{'N/A':>10}"
        s_str = f"{mode.entropy:12.6f}" if mode.entropy is not None else f"{'N/A':>12}"
        marker = " *" if top is not None and mode.mode_id == top.mode_id else ""
        print(
            f"{mode.mode_id:>5}  {mode.rank:>4}  {mode.n_poses:>5}  "
            f"{cf_str}  {fe_str}  {h_str}  {s_str}{marker}"
        )

    if top is not None:
        print(f"\n* = top-ranked mode (mode_id={top.mode_id})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
