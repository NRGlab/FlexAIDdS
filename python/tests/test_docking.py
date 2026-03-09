"""Tests for python/flexaidds/docking.py.

Config-file parsing and pure-Python binding mode / population APIs are tested
here without requiring the compiled C++ extension.
"""

from __future__ import annotations

import pytest


# ─────────────────────────────────────────────────────────────────────────────
# Docking – config file parsing
# ─────────────────────────────────────────────────────────────────────────────

def test_docking_missing_config(tmp_path):
    """Docking raises FileNotFoundError for a non-existent config."""
    from flexaidds.docking import Docking

    with pytest.raises(FileNotFoundError):
        Docking(str(tmp_path / "nonexistent.inp"))


def test_docking_parse_string_keys(flexaid_config_file):
    """String-valued keywords (PDBNAM, INPLIG, METOPT) are parsed correctly."""
    from flexaidds.docking import Docking

    d = Docking(str(flexaid_config_file))
    assert d.receptor == "/data/receptor.pdb"
    assert d.ligand   == "/data/ligand.mol2"
    assert d.optimization_method == "GA"


def test_docking_parse_int_keys(flexaid_config_file):
    """Integer keywords (TEMPER, NRGOUT) are stored as int."""
    from flexaidds.docking import Docking

    d = Docking(str(flexaid_config_file))
    assert d.temperature == 300
    assert isinstance(d.temperature, int)


def test_docking_parse_list_keys(flexaid_config_file):
    """Keywords that can repeat (OPTIMZ) are stored as a list."""
    from flexaidds.docking import Docking

    d = Docking(str(flexaid_config_file))
    assert d._config["OPTIMZ"] == ["-12.0", "-11.5"]


def test_docking_parse_flag_keys(flexaid_config_file):
    """Flag-only keywords (EXCHET) are stored as True."""
    from flexaidds.docking import Docking

    d = Docking(str(flexaid_config_file))
    assert d._config.get("EXCHET") is True


def test_docking_parse_float_keys(flexaid_config_file):
    """Float keywords (ACSWEI) are stored as float."""
    from flexaidds.docking import Docking

    d = Docking(str(flexaid_config_file))
    assert isinstance(d._config["ACSWEI"], float)
    assert abs(d._config["ACSWEI"] - 0.75) < 1e-9


def test_docking_run_not_implemented(flexaid_config_file):
    """Docking.run() raises NotImplementedError (full GA not yet wired)."""
    from flexaidds.docking import Docking

    d = Docking(str(flexaid_config_file))
    with pytest.raises(NotImplementedError):
        d.run()


def test_docking_repr(flexaid_config_file):
    """Docking.__repr__ includes the config filename."""
    from flexaidds.docking import Docking

    d = Docking(str(flexaid_config_file))
    assert "test.inp" in repr(d)


# ─────────────────────────────────────────────────────────────────────────────
# BindingMode – pure Python path
# ─────────────────────────────────────────────────────────────────────────────

def test_binding_mode_no_cpp():
    """BindingMode with no C++ object returns sentinel / zero values."""
    from flexaidds.docking import BindingMode

    mode = BindingMode(cpp_binding_mode=None)
    assert mode.n_poses == 0
    assert mode.free_energy == float("inf")
    assert mode.entropy == 0.0
    assert len(mode) == 0


# ─────────────────────────────────────────────────────────────────────────────
# BindingPopulation – pure Python path
# ─────────────────────────────────────────────────────────────────────────────

def test_binding_population_empty():
    """BindingPopulation starts empty."""
    from flexaidds.docking import BindingPopulation

    pop = BindingPopulation()
    assert pop.n_modes == 0
    assert len(pop) == 0


def test_binding_population_rank_by_free_energy():
    """rank_by_free_energy returns modes in ascending F order."""
    from flexaidds.docking import BindingMode, BindingPopulation

    class _FakeMode(BindingMode):
        def __init__(self, F):
            super().__init__(cpp_binding_mode=None)
            self._F = F
        @property
        def free_energy(self):
            return self._F

    pop = BindingPopulation()
    for F in [0.5, -1.2, -3.0, -0.1]:
        pop.add_mode(_FakeMode(F))

    ranked = pop.rank_by_free_energy()
    Fs = [m.free_energy for m in ranked]
    assert Fs == sorted(Fs)


def test_binding_population_iter_and_getitem():
    """BindingPopulation supports iteration and indexing."""
    from flexaidds.docking import BindingMode, BindingPopulation

    pop = BindingPopulation()
    modes = [BindingMode() for _ in range(3)]
    for m in modes:
        pop.add_mode(m)

    assert len(list(pop)) == 3
    assert pop[0] is modes[0]
    assert pop[2] is modes[2]
