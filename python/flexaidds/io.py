"""File I/O utilities for FlexAID∆S.

Parsers and writers for the file formats consumed and produced by the
FlexAID docking engine:

    - FlexAID .inp configuration files  (read)
    - RRD result files                  (read)
    - PDB coordinate files              (read / write)
    - MOL2 ligand files                 (read)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple


# ─────────────────────────────────────────────────────────────────────────────
# PDB atoms
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Atom:
    """A single ATOM / HETATM record from a PDB file.

    Attributes:
        serial:    PDB serial number (1-based).
        name:      Atom name (e.g. ``"CA"``).
        alt_loc:   Alternate location indicator (usually ``" "``).
        res_name:  Residue name (e.g. ``"ALA"``).
        chain_id:  Chain identifier.
        res_seq:   Residue sequence number.
        i_code:    Insertion code.
        x, y, z:  Cartesian coordinates (Å).
        occupancy: Occupancy (0.0–1.0).
        b_factor:  Temperature factor.
        element:   Element symbol.
        hetatm:    True for HETATM records.
    """
    serial: int = 0
    name: str = ""
    alt_loc: str = " "
    res_name: str = ""
    chain_id: str = " "
    res_seq: int = 0
    i_code: str = " "
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    occupancy: float = 1.0
    b_factor: float = 0.0
    element: str = ""
    hetatm: bool = False

    @property
    def coords(self) -> Tuple[float, float, float]:
        """(x, y, z) tuple."""
        return (self.x, self.y, self.z)

    def to_pdb_line(self) -> str:
        """Format as a fixed-width PDB ATOM/HETATM line."""
        record = "HETATM" if self.hetatm else "ATOM  "
        return (
            f"{record}{self.serial:5d} {self.name:<4s}{self.alt_loc}"
            f"{self.res_name:<3s} {self.chain_id}{self.res_seq:4d}{self.i_code}   "
            f"{self.x:8.3f}{self.y:8.3f}{self.z:8.3f}"
            f"{self.occupancy:6.2f}{self.b_factor:6.2f}          "
            f"{self.element:>2s}\n"
        )

    def __repr__(self) -> str:
        return (
            f"<Atom {self.serial} {self.name!r} "
            f"{self.res_name}{self.res_seq} ({self.x:.2f},{self.y:.2f},{self.z:.2f})>"
        )


def read_pdb(path: str | Path) -> List[Atom]:
    """Parse ATOM and HETATM records from a PDB file.

    Args:
        path: Path to the PDB file.

    Returns:
        List of :class:`Atom` objects in file order.
    """
    atoms: List[Atom] = []
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                atom = Atom(
                    serial   = int(line[6:11]),
                    name     = line[12:16].strip(),
                    alt_loc  = line[16],
                    res_name = line[17:20].strip(),
                    chain_id = line[21],
                    res_seq  = int(line[22:26]),
                    i_code   = line[26],
                    x        = float(line[30:38]),
                    y        = float(line[38:46]),
                    z        = float(line[46:54]),
                    occupancy= float(line[54:60]) if len(line) > 54 else 1.0,
                    b_factor = float(line[60:66]) if len(line) > 60 else 0.0,
                    element  = line[76:78].strip() if len(line) > 76 else "",
                    hetatm   = (rec == "HETATM"),
                )
            except (ValueError, IndexError):
                continue
            atoms.append(atom)
    return atoms


def write_pdb(atoms: List[Atom], path: str | Path) -> None:
    """Write a list of atoms to a PDB file.

    Args:
        atoms: Atoms to write.
        path:  Output path.
    """
    with open(path, "w") as fh:
        for atom in atoms:
            fh.write(atom.to_pdb_line())
        fh.write("END\n")


# ─────────────────────────────────────────────────────────────────────────────
# RRD result files
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class RRDPose:
    """A single docked pose from a FlexAID RRD result file.

    Attributes:
        pose_id:    Pose number (1-based).
        cf_score:   Contact function score (kcal/mol; lower is better).
        rmsd:       RMSD to reference (Å), if present in the file.
        mode_id:    Binding mode cluster index.
        chromosome: Raw GA chromosome string.
    """
    pose_id: int = 0
    cf_score: float = 0.0
    rmsd: Optional[float] = None
    mode_id: int = 0
    chromosome: str = ""

    def __repr__(self) -> str:
        rmsd_str = f" rmsd={self.rmsd:.2f}" if self.rmsd is not None else ""
        return (
            f"<RRDPose {self.pose_id} mode={self.mode_id} "
            f"CF={self.cf_score:.4f}{rmsd_str}>"
        )


def read_rrd(path: str | Path) -> List[RRDPose]:
    """Parse a FlexAID RRD (ranked results data) file.

    The RRD format is whitespace-delimited:
        pose_id  cf_score  [rmsd]  mode_id  chromosome

    Lines starting with ``#`` are ignored.

    Args:
        path: Path to the ``.rrd`` file.

    Returns:
        List of :class:`RRDPose` objects.
    """
    poses: List[RRDPose] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                pose_id  = int(parts[0])
                cf_score = float(parts[1])
                # Detect RMSD column: try to parse parts[2] as int.
                # If it succeeds it's the mode_id (no RMSD column);
                # if it fails (decimal point) it's a float RMSD.
                try:
                    mode_id = int(parts[2])
                    rmsd    = None
                    chrom   = " ".join(parts[3:])
                except ValueError:
                    rmsd    = float(parts[2])
                    mode_id = int(parts[3])
                    chrom   = " ".join(parts[4:])
            except (ValueError, IndexError):
                continue
            poses.append(RRDPose(
                pose_id=pose_id,
                cf_score=cf_score,
                rmsd=rmsd,
                mode_id=mode_id,
                chromosome=chrom,
            ))
    return poses


# ─────────────────────────────────────────────────────────────────────────────
# MOL2 ligand files (minimal ATOM block reader)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Mol2Atom:
    """A single atom from a Tripos MOL2 @<TRIPOS>ATOM block.

    Attributes:
        atom_id:   1-based integer atom ID.
        atom_name: Atom name string.
        x, y, z:  Cartesian coordinates (Å).
        atom_type: Tripos atom type (e.g. ``"C.3"``).
        res_id:    Residue/substructure ID.
        res_name:  Residue/substructure name.
        charge:    Partial charge (e).
    """
    atom_id: int = 0
    atom_name: str = ""
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    atom_type: str = ""
    res_id: int = 1
    res_name: str = "LIG"
    charge: float = 0.0

    @property
    def element(self) -> str:
        """Element symbol derived from Tripos atom type (e.g. ``"C.3"`` → ``"C"``)."""
        return self.atom_type.split(".")[0]

    @property
    def coords(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)

    def __repr__(self) -> str:
        return f"<Mol2Atom {self.atom_id} {self.atom_name!r} ({self.x:.2f},{self.y:.2f},{self.z:.2f})>"


def read_mol2(path: str | Path) -> List[Mol2Atom]:
    """Parse the @<TRIPOS>ATOM section from a MOL2 file.

    Args:
        path: Path to the ``.mol2`` file.

    Returns:
        List of :class:`Mol2Atom` objects.
    """
    atoms: List[Mol2Atom] = []
    in_atom_block = False
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("@<TRIPOS>ATOM"):
                in_atom_block = True
                continue
            if stripped.startswith("@<TRIPOS>"):
                in_atom_block = False
                continue
            if not in_atom_block or not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 6:
                continue
            try:
                atoms.append(Mol2Atom(
                    atom_id   = int(parts[0]),
                    atom_name = parts[1],
                    x         = float(parts[2]),
                    y         = float(parts[3]),
                    z         = float(parts[4]),
                    atom_type = parts[5],
                    res_id    = int(parts[6])   if len(parts) > 6 else 1,
                    res_name  = parts[7]        if len(parts) > 7 else "LIG",
                    charge    = float(parts[8]) if len(parts) > 8 else 0.0,
                ))
            except (ValueError, IndexError):
                continue
    return atoms
