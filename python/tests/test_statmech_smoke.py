import pytest
import flexaidds as fd

_core_available: bool
try:
    fd.StatMechEngine(300.0)
    _core_available = True
except Exception:
    _core_available = False


@pytest.mark.skipif(not _core_available, reason="C++ _core extension not built")
def test_statmech_smoke():
    engine = fd.StatMechEngine(300.0)
    engine.add_sample(-7.0)
    engine.add_sample(-6.0)
    thermo = engine.compute()

    assert hasattr(thermo, "free_energy")
    assert hasattr(thermo, "mean_energy")
    assert hasattr(thermo, "entropy")
    assert hasattr(thermo, "heat_capacity")
    assert hasattr(thermo, "std_energy")
