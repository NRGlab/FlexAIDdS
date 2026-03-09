"""Tests for python/flexaidds/io.py.

All tests are pure-Python (no C++ extension required).
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import pytest


# ─────────────────────────────────────────────────────────────────────────────
# PDB I/O
# ─────────────────────────────────────────────────────────────────────────────

def test_read_pdb_atom_count(simple_pdb_file):
    """read_pdb returns all ATOM and HETATM records."""
    from flexaidds.io import read_pdb

    atoms = read_pdb(simple_pdb_file)
    assert len(atoms) == 3


def test_read_pdb_coords(simple_pdb_file):
    """Coordinates are parsed to float with correct values."""
    from flexaidds.io import read_pdb

    atoms = read_pdb(simple_pdb_file)
    assert math.isclose(atoms[0].x, 1.0, abs_tol=1e-4)
    assert math.isclose(atoms[0].y, 2.0, abs_tol=1e-4)
    assert math.isclose(atoms[0].z, 3.0, abs_tol=1e-4)


def test_read_pdb_hetatm_flag(simple_pdb_file):
    """HETATM records set hetatm=True."""
    from flexaidds.io import read_pdb

    atoms = read_pdb(simple_pdb_file)
    assert atoms[0].hetatm is False
    assert atoms[2].hetatm is True


def test_read_pdb_residue_fields(simple_pdb_file):
    """Residue name and chain ID are parsed correctly."""
    from flexaidds.io import read_pdb

    atoms = read_pdb(simple_pdb_file)
    assert atoms[0].res_name == "ALA"
    assert atoms[0].chain_id == "A"
    assert atoms[2].chain_id == "B"


def test_atom_coords_property(simple_pdb_file):
    """Atom.coords returns a 3-tuple."""
    from flexaidds.io import read_pdb

    atom = read_pdb(simple_pdb_file)[1]
    x, y, z = atom.coords
    assert math.isclose(x, 2.0, abs_tol=1e-4)


def test_write_pdb_roundtrip(simple_pdb_file, tmp_path):
    """write_pdb + read_pdb roundtrip preserves atom count and coordinates."""
    from flexaidds.io import read_pdb, write_pdb

    atoms = read_pdb(simple_pdb_file)
    out = tmp_path / "out.pdb"
    write_pdb(atoms, out)

    atoms2 = read_pdb(out)
    assert len(atoms2) == len(atoms)
    for a, b in zip(atoms, atoms2):
        assert math.isclose(a.x, b.x, abs_tol=1e-3)
        assert math.isclose(a.y, b.y, abs_tol=1e-3)
        assert math.isclose(a.z, b.z, abs_tol=1e-3)


# ─────────────────────────────────────────────────────────────────────────────
# RRD I/O
# ─────────────────────────────────────────────────────────────────────────────

def test_read_rrd_pose_count(simple_rrd_file):
    """read_rrd returns correct number of poses."""
    from flexaidds.io import read_rrd

    poses = read_rrd(simple_rrd_file)
    assert len(poses) == 3


def test_read_rrd_cf_scores(simple_rrd_file):
    """CF scores are parsed as floats."""
    from flexaidds.io import read_rrd

    poses = read_rrd(simple_rrd_file)
    assert math.isclose(poses[0].cf_score, -12.5, abs_tol=1e-9)
    assert math.isclose(poses[2].cf_score, -10.2, abs_tol=1e-9)


def test_read_rrd_mode_ids(simple_rrd_file):
    """Mode IDs are parsed correctly."""
    from flexaidds.io import read_rrd

    poses = read_rrd(simple_rrd_file)
    assert poses[0].mode_id == 1
    assert poses[2].mode_id == 2


def test_read_rrd_no_rmsd(simple_rrd_file):
    """RRD without RMSD column sets rmsd=None."""
    from flexaidds.io import read_rrd

    poses = read_rrd(simple_rrd_file)
    assert poses[0].rmsd is None


def test_read_rrd_with_rmsd(rrd_with_rmsd_file):
    """RRD with RMSD column parses rmsd correctly."""
    from flexaidds.io import read_rrd

    poses = read_rrd(rrd_with_rmsd_file)
    assert poses[0].rmsd is not None
    assert math.isclose(poses[0].rmsd, 0.45, abs_tol=1e-9)
    assert math.isclose(poses[1].rmsd, 1.20, abs_tol=1e-9)


def test_read_rrd_repr(simple_rrd_file):
    """RRDPose repr contains pose ID and mode ID."""
    from flexaidds.io import read_rrd

    pose = read_rrd(simple_rrd_file)[0]
    r = repr(pose)
    assert "RRDPose" in r
    assert "mode=1" in r


# ─────────────────────────────────────────────────────────────────────────────
# MOL2 I/O
# ─────────────────────────────────────────────────────────────────────────────

def test_read_mol2_atom_count(simple_mol2_file):
    """read_mol2 parses both atoms from the ATOM block."""
    from flexaidds.io import read_mol2

    atoms = read_mol2(simple_mol2_file)
    assert len(atoms) == 2


def test_read_mol2_coords(simple_mol2_file):
    """MOL2 atom coordinates are parsed correctly."""
    from flexaidds.io import read_mol2

    atoms = read_mol2(simple_mol2_file)
    assert math.isclose(atoms[0].x, 0.0, abs_tol=1e-4)
    assert math.isclose(atoms[1].x, 1.4, abs_tol=1e-4)


def test_read_mol2_atom_type(simple_mol2_file):
    """Tripos atom types are preserved."""
    from flexaidds.io import read_mol2

    atoms = read_mol2(simple_mol2_file)
    assert atoms[0].atom_type == "C.3"
    assert atoms[1].atom_type == "O.3"


def test_read_mol2_element(simple_mol2_file):
    """Element is derived from Tripos atom type."""
    from flexaidds.io import read_mol2

    atoms = read_mol2(simple_mol2_file)
    assert atoms[0].element == "C"
    assert atoms[1].element == "O"


def test_read_mol2_charge(simple_mol2_file):
    """Partial charges are parsed as float."""
    from flexaidds.io import read_mol2

    atoms = read_mol2(simple_mol2_file)
    assert math.isclose(atoms[1].charge, -0.3982, abs_tol=1e-4)
