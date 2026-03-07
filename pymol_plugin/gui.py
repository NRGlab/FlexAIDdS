"""FlexAID∆S PyMOL GUI Panel.

Qt-based interface for docking result visualization and analysis.
"""

import os
from pathlib import Path

try:
    from pymol.Qt import QtWidgets, QtCore
    from pymol import cmd
except ImportError:
    raise ImportError("PyMOL Qt bindings not available")

from .visualization import (
    load_binding_modes,
    show_pose_ensemble,
    color_by_boltzmann_weight,
    show_thermodynamics,
    export_to_nrgsuite,
    _loaded_modes,
    _temperature_K,
)


class FlexAIDSPanel(QtWidgets.QDialog):
    """Main GUI panel for FlexAID∆S plugin."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("FlexAID∆S: Entropy-Driven Docking")
        self.setMinimumWidth(400)
        self.setMinimumHeight(500)

        self._setup_ui()
        self._connect_signals()

        # Data state — mode_name strings parallel to QListWidget rows
        self._mode_names: list = []

    def _setup_ui(self):
        """Construct GUI layout."""
        layout = QtWidgets.QVBoxLayout(self)

        # ─── File loading section ───
        file_group = QtWidgets.QGroupBox("Load Docking Results")
        file_layout = QtWidgets.QHBoxLayout()

        self.file_path_edit = QtWidgets.QLineEdit()
        self.file_path_edit.setPlaceholderText("Select FlexAID output directory...")

        self.browse_btn = QtWidgets.QPushButton("Browse")
        self.load_btn = QtWidgets.QPushButton("Load")
        self.load_btn.setEnabled(False)

        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.browse_btn)
        file_layout.addWidget(self.load_btn)
        file_group.setLayout(file_layout)
        layout.addWidget(file_group)

        # ─── Binding mode list ───
        mode_group = QtWidgets.QGroupBox("Binding Modes")
        mode_layout = QtWidgets.QVBoxLayout()

        self.mode_list = QtWidgets.QListWidget()
        self.mode_list.setSelectionMode(QtWidgets.QListWidget.SingleSelection)
        mode_layout.addWidget(self.mode_list)

        mode_group.setLayout(mode_layout)
        layout.addWidget(mode_group)

        # ─── Thermodynamic properties display ───
        thermo_group = QtWidgets.QGroupBox("Thermodynamics")
        thermo_layout = QtWidgets.QFormLayout()

        self.free_energy_label = QtWidgets.QLabel("-")
        self.enthalpy_label = QtWidgets.QLabel("-")
        self.entropy_label = QtWidgets.QLabel("-")
        self.entropy_term_label = QtWidgets.QLabel("-")
        self.n_poses_label = QtWidgets.QLabel("-")

        thermo_layout.addRow("ΔG (kcal/mol):", self.free_energy_label)
        thermo_layout.addRow("ΔH (kcal/mol):", self.enthalpy_label)
        thermo_layout.addRow("S (kcal/mol·K):", self.entropy_label)
        thermo_layout.addRow("TΔS (kcal/mol):", self.entropy_term_label)
        thermo_layout.addRow("# Poses:", self.n_poses_label)

        thermo_group.setLayout(thermo_layout)
        layout.addWidget(thermo_group)

        # ─── Visualization controls ───
        viz_group = QtWidgets.QGroupBox("Visualization")
        viz_layout = QtWidgets.QVBoxLayout()

        self.show_ensemble_btn = QtWidgets.QPushButton("Show Pose Ensemble")
        self.show_ensemble_btn.setEnabled(False)

        self.color_boltzmann_btn = QtWidgets.QPushButton("Color by Boltzmann Weight")
        self.color_boltzmann_btn.setEnabled(False)

        self.show_representative_btn = QtWidgets.QPushButton("Show Representative Only")
        self.show_representative_btn.setEnabled(False)

        viz_layout.addWidget(self.show_ensemble_btn)
        viz_layout.addWidget(self.color_boltzmann_btn)
        viz_layout.addWidget(self.show_representative_btn)

        viz_group.setLayout(viz_layout)
        layout.addWidget(viz_group)

        # ─── NRGSuite integration ───
        nrg_group = QtWidgets.QGroupBox("NRGSuite Integration")
        nrg_layout = QtWidgets.QVBoxLayout()

        self.launch_nrgsuite_btn = QtWidgets.QPushButton("Launch NRGSuite")
        self.export_to_nrg_btn = QtWidgets.QPushButton("Export to NRGSuite Format")
        self.export_to_nrg_btn.setEnabled(False)

        nrg_layout.addWidget(self.launch_nrgsuite_btn)
        nrg_layout.addWidget(self.export_to_nrg_btn)

        nrg_group.setLayout(nrg_layout)
        layout.addWidget(nrg_group)

        # ─── Close button ───
        close_btn = QtWidgets.QPushButton("Close")
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)

    def _connect_signals(self):
        """Wire up button click handlers."""
        self.browse_btn.clicked.connect(self._browse_directory)
        self.file_path_edit.textChanged.connect(
            lambda t: self.load_btn.setEnabled(bool(t.strip()))
        )
        self.load_btn.clicked.connect(self._load_results)
        self.mode_list.itemSelectionChanged.connect(self._on_mode_selected)
        self.show_ensemble_btn.clicked.connect(self._show_pose_ensemble)
        self.color_boltzmann_btn.clicked.connect(self._color_by_boltzmann)
        self.show_representative_btn.clicked.connect(self._show_representative)
        self.launch_nrgsuite_btn.clicked.connect(self._launch_nrgsuite)
        self.export_to_nrg_btn.clicked.connect(self._export_to_nrgsuite)

    # ── slots ────────────────────────────────────────────────────────────────

    def _browse_directory(self):
        """Open file dialog to select FlexAID output directory."""
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Select FlexAID Output Directory"
        )
        if directory:
            self.file_path_edit.setText(directory)

    def _load_results(self):
        """Load docking results from the selected output directory."""
        output_dir = Path(self.file_path_edit.text().strip())
        if not output_dir.exists():
            QtWidgets.QMessageBox.warning(
                self, "Error", f"Directory not found: {output_dir}"
            )
            return

        # Delegate to visualization module (parses PDBs, .cad, computes thermo)
        load_binding_modes(str(output_dir), temperature=_temperature_K)

        if not _loaded_modes:
            QtWidgets.QMessageBox.warning(
                self, "No Results",
                f"No FlexAID result files found in:\n{output_dir}\n\n"
                "Expected: result_*.pdb files."
            )
            return

        # Populate list widget with ranked modes (sort by free energy)
        self.mode_list.clear()
        self._mode_names.clear()

        sorted_modes = sorted(
            _loaded_modes.items(),
            key=lambda kv: (
                kv[1].free_energy if kv[1].free_energy is not None else float("inf")
            ),
        )

        for mode_name, rec in sorted_modes:
            dg_str = (
                f"{rec.free_energy:.2f} kcal/mol"
                if rec.free_energy is not None
                else "N/A"
            )
            label = f"{mode_name}: ΔG = {dg_str}  ({rec.frequency} poses)"
            self.mode_list.addItem(label)
            self._mode_names.append(mode_name)

        # Enable controls
        for btn in (
            self.show_ensemble_btn,
            self.color_boltzmann_btn,
            self.show_representative_btn,
            self.export_to_nrg_btn,
        ):
            btn.setEnabled(True)

        # Auto-select first mode
        if self.mode_list.count():
            self.mode_list.setCurrentRow(0)

        QtWidgets.QMessageBox.information(
            self,
            "Loaded",
            f"Loaded {len(_loaded_modes)} binding modes from {output_dir.name}.",
        )

    def _on_mode_selected(self):
        """Update thermodynamics panel when a binding mode is selected."""
        row = self.mode_list.currentRow()
        if row < 0 or row >= len(self._mode_names):
            return

        mode_name = self._mode_names[row]
        rec = _loaded_modes.get(mode_name)
        if rec is None:
            return

        T = _temperature_K
        entropy_term = (rec.entropy * T) if rec.entropy is not None else None

        def _fmt(val, fmt=".4f"):
            return f"{val:{fmt}}" if val is not None else "N/A"

        self.free_energy_label.setText(_fmt(rec.free_energy))
        self.enthalpy_label.setText(_fmt(rec.enthalpy))
        self.entropy_label.setText(_fmt(rec.entropy, ".6f"))
        self.entropy_term_label.setText(_fmt(entropy_term))
        self.n_poses_label.setText(str(rec.frequency))

    def _selected_mode_name(self) -> str | None:
        """Return the currently selected mode name, or None."""
        row = self.mode_list.currentRow()
        if row < 0 or row >= len(self._mode_names):
            QtWidgets.QMessageBox.warning(
                self, "No Mode Selected", "Please select a binding mode first."
            )
            return None
        return self._mode_names[row]

    def _show_pose_ensemble(self):
        """Render all poses in the selected binding mode."""
        mode_name = self._selected_mode_name()
        if mode_name:
            show_pose_ensemble(mode_name, show_all=True)

    def _color_by_boltzmann(self):
        """Color poses by Boltzmann weight (blue = low, red = high)."""
        mode_name = self._selected_mode_name()
        if mode_name:
            color_by_boltzmann_weight(mode_name)

    def _show_representative(self):
        """Show only the representative pose (highest Boltzmann weight)."""
        mode_name = self._selected_mode_name()
        if mode_name:
            show_pose_ensemble(mode_name, show_all=False)

    def _launch_nrgsuite(self):
        """Launch NRGSuite plugin (if installed)."""
        try:
            cmd.do("nrgsuite")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(
                self,
                "NRGSuite Not Available",
                f"Could not launch NRGSuite: {exc}\n\n"
                "Install from: https://github.com/NRGlab/NRGsuite",
            )

    def _export_to_nrgsuite(self):
        """Export binding modes to NRGSuite-compatible TSV format."""
        output_dir = self.file_path_edit.text().strip()
        if not output_dir:
            QtWidgets.QMessageBox.warning(
                self, "No Directory", "Load a results directory first."
            )
            return

        nrg_file, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save NRGSuite Export",
            str(Path(output_dir) / "nrgsuite_export.txt"),
            "Text Files (*.txt);;All Files (*)",
        )
        if not nrg_file:
            return

        export_to_nrgsuite(output_dir, nrg_file)
        QtWidgets.QMessageBox.information(
            self, "Exported", f"NRGSuite file written to:\n{nrg_file}"
        )
