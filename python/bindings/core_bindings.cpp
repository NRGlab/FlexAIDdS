// core_bindings.cpp — pybind11 bindings for FlexAID∆S C++ core
//
// Exposes:
//   - statmech::StatMechEngine
//   - statmech::Thermodynamics
//   - statmech::BoltzmannLUT
//   - BindingMode (basic interface)
//   - BindingPopulation (basic interface)
//
// Build: See python/setup.py and CMakeLists.txt with -DBUILD_PYTHON_BINDINGS=ON

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "../../LIB/statmech.h"
#include "../../LIB/BindingMode.h"

namespace py = pybind11;
using namespace statmech;

// ──────────────────────────────────────────────────────────────────────────────
// Helper: Convert C++ std::vector to NumPy array (zero-copy view when possible)
// ──────────────────────────────────────────────────────────────────────────────

template <typename T>
py::array_t<T> to_numpy(const std::vector<T>& vec) {
    return py::array_t<T>(vec.size(), vec.data());
}

// ──────────────────────────────────────────────────────────────────────────────
// Module definition
// ──────────────────────────────────────────────────────────────────────────────

PYBIND11_MODULE(_core, m) {
    m.doc() = "FlexAID∆S C++ core bindings: statistical mechanics and docking engine";
    
    // Physical constants
    m.attr("kB_kcal") = kB_kcal;
    m.attr("kB_SI")   = kB_SI;
    
    // ═══════════════════════════════════════════════════════════════════════
    // Thermodynamics data structures
    // ═══════════════════════════════════════════════════════════════════════
    
    py::class_<State>(m, "State", "Micr ostate with energy and degeneracy")
        .def(py::init<>())
        .def_readwrite("energy", &State::energy, "Energy in kcal/mol")
        .def_readwrite("count",  &State::count,  "Degeneracy/multiplicity")
        .def("__repr__", [](const State& s) {
            return "<State energy=" + std::to_string(s.energy) + 
                   " count=" + std::to_string(s.count) + ">";
        });
    
    py::class_<Thermodynamics>(m, "Thermodynamics", 
        "Complete thermodynamic properties of an ensemble")
        .def(py::init<>())
        .def_readwrite("temperature",    &Thermodynamics::temperature,    "K")
        .def_readwrite("log_Z",          &Thermodynamics::log_Z,          "ln(partition function)")
        .def_readwrite("free_energy",    &Thermodynamics::free_energy,    "Helmholtz F (kcal/mol)")
        .def_readwrite("mean_energy",    &Thermodynamics::mean_energy,    "⟨E⟩ (kcal/mol)")
        .def_readwrite("mean_energy_sq", &Thermodynamics::mean_energy_sq, "⟨E²⟩")
        .def_readwrite("heat_capacity",  &Thermodynamics::heat_capacity,  "Cv (kcal mol⁻¹ K⁻²)")
        .def_readwrite("entropy",        &Thermodynamics::entropy,        "S (kcal mol⁻¹ K⁻¹)")
        .def_readwrite("std_energy",     &Thermodynamics::std_energy,     "σ_E (kcal/mol)")
        .def("__repr__", [](const Thermodynamics& t) {
            char buf[256];
            snprintf(buf, sizeof(buf),
                "<Thermodynamics T=%.1fK F=%.3f H=%.3f S=%.5f Cv=%.3f>",
                t.temperature, t.free_energy, t.mean_energy, 
                t.entropy, t.heat_capacity);
            return std::string(buf);
        });
    
    py::class_<Replica>(m, "Replica", "Parallel tempering replica")
        .def(py::init<>())
        .def_readwrite("id",              &Replica::id)
        .def_readwrite("temperature",     &Replica::temperature)
        .def_readwrite("beta",            &Replica::beta)
        .def_readwrite("current_energy",  &Replica::current_energy);
    
    py::class_<WHAMBin>(m, "WHAMBin", "WHAM histogram bin with free energy")
        .def(py::init<>())
        .def_readwrite("coord_center",  &WHAMBin::coord_center)
        .def_readwrite("count",         &WHAMBin::count)
        .def_readwrite("free_energy",   &WHAMBin::free_energy);
    
    py::class_<TIPoint>(m, "TIPoint", "Thermodynamic integration data point")
        .def(py::init<>())
        .def_readwrite("lambda",       &TIPoint::lambda)
        .def_readwrite("dV_dlambda",   &TIPoint::dV_dlambda);
    
    // ═══════════════════════════════════════════════════════════════════════
    // StatMechEngine: core thermodynamics calculator
    // ═══════════════════════════════════════════════════════════════════════
    
    py::class_<StatMechEngine>(m, "StatMechEngine", 
        "Statistical mechanics engine for conformational ensembles")
        .def(py::init<double>(), 
            py::arg("temperature_K") = 300.0,
            "Initialize engine at given temperature (default 300K)")
        
        // ─── Ensemble construction ───
        .def("add_sample", &StatMechEngine::add_sample,
            py::arg("energy"), py::arg("multiplicity") = 1,
            "Add a sampled configuration with energy (kcal/mol) and multiplicity")
        .def("clear", &StatMechEngine::clear, "Remove all samples")
        
        // ─── Thermodynamic analysis ───
        .def("compute", &StatMechEngine::compute,
            "Compute full thermodynamics (F, S, H, Cv, etc.)")
        .def("boltzmann_weights", &StatMechEngine::boltzmann_weights,
            "Return Boltzmann weights for all samples (same order as insertion)")
        .def("delta_G", &StatMechEngine::delta_G,
            py::arg("reference"),
            "Compute ΔG relative to another ensemble")
        
        // ─── Advanced methods ───
        .def_static("init_replicas", &StatMechEngine::init_replicas,
            py::arg("temperatures"),
            "Initialize parallel tempering replicas at given temperatures")
        .def_static("attempt_swap", 
            [](Replica& a, Replica& b) {
                // Python RNG not compatible with std::mt19937, use C++ RNG
                static std::mt19937 rng{std::random_device{}()};
                return StatMechEngine::attempt_swap(a, b, rng);
            },
            py::arg("replica_a"), py::arg("replica_b"),
            "Attempt Metropolis swap between two replicas (returns True if accepted)")
        .def_static("wham", &StatMechEngine::wham,
            py::arg("energies"), py::arg("coordinates"), 
            py::arg("temperature"), py::arg("n_bins"),
            py::arg("max_iter") = 1000, py::arg("tolerance") = 1e-6,
            "WHAM free energy profile along a reaction coordinate")
        .def_static("thermodynamic_integration", 
            &StatMechEngine::thermodynamic_integration,
            py::arg("points"),
            "Compute ΔG via thermodynamic integration (trapezoidal rule)")
        .def_static("helmholtz", &StatMechEngine::helmholtz,
            py::arg("energies"), py::arg("temperature"),
            "Compute Helmholtz free energy from raw energy vector")
        
        // ─── Properties ───
        .def_property_readonly("temperature", &StatMechEngine::temperature)
        .def_property_readonly("beta", &StatMechEngine::beta)
        .def_property_readonly("size", &StatMechEngine::size)
        .def("__len__", &StatMechEngine::size)
        .def("__repr__", [](const StatMechEngine& e) {
            return "<StatMechEngine T=" + std::to_string(e.temperature()) + 
                   "K n_samples=" + std::to_string(e.size()) + ">";
        });
    
    // ═══════════════════════════════════════════════════════════════════════
    // BoltzmannLUT: fast lookup table
    // ═══════════════════════════════════════════════════════════════════════
    
    py::class_<BoltzmannLUT>(m, "BoltzmannLUT",
        "Pre-tabulated Boltzmann factors for O(1) inner-loop evaluation")
        .def(py::init<double, double, double, int>(),
            py::arg("beta"), py::arg("e_min"), py::arg("e_max"),
            py::arg("n_bins") = 10000,
            "Initialize LUT for energy range [e_min, e_max]")
        .def("__call__", &BoltzmannLUT::operator(),
            py::arg("energy"),
            "Look up exp(-β E) for given energy");
    
    // ═══════════════════════════════════════════════════════════════════════
    // BindingMode: pose cluster with thermodynamic scoring
    // ═══════════════════════════════════════════════════════════════════════
    
    py::class_<BindingMode>(m, "BindingMode",
        "Binding mode: cluster of docked poses with thermodynamic analysis")
        // Legacy interface (backward compatibility)
        .def("compute_energy", &BindingMode::compute_energy,
            "Helmholtz free energy F = H - TS (kcal/mol)")
        .def("compute_entropy", &BindingMode::compute_entropy,
            "Configurational entropy S (kcal mol⁻¹ K⁻¹)")
        .def("compute_enthalpy", &BindingMode::compute_enthalpy,
            "Boltzmann-weighted average energy ⟨E⟩ (kcal/mol)")
        
        // New StatMech API
        .def("get_thermodynamics", &BindingMode::get_thermodynamics,
            "Full thermodynamic properties (F, S, H, Cv, σ_E)")
        .def("get_free_energy", &BindingMode::get_free_energy,
            "Alias for compute_energy()")
        .def("get_heat_capacity", &BindingMode::get_heat_capacity,
            "Heat capacity Cv (kcal mol⁻¹ K⁻²)")
        .def("get_boltzmann_weights", &BindingMode::get_boltzmann_weights,
            "Boltzmann weights for all poses in this mode")
        .def("get_BindingMode_size", &BindingMode::get_BindingMode_size,
            "Number of poses in this binding mode")
        .def("__len__", &BindingMode::get_BindingMode_size)
        .def("__repr__", [](const BindingMode& m) {
            auto thermo = m.get_thermodynamics();
            char buf[256];
            snprintf(buf, sizeof(buf),
                "<BindingMode n_poses=%d F=%.3f H=%.3f S=%.5f>",
                m.get_BindingMode_size(), 
                thermo.free_energy, thermo.mean_energy, thermo.entropy);
            return std::string(buf);
        });
    
    // Note: BindingPopulation requires GA/docking infrastructure.
    // Full integration deferred to Phase 2 (see python/flexaidds/docking.py
    // for high-level wrapper when C++ FlexAID docking engine is exposed).
}
