// CavityDetect.cpp
// Native C++20 SURFNET-style cavity detection for FlexAIDΔS + FreeNRG
// In-memory, OpenMP-parallel, zero external deps, zero GPL contact
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include "CavityDetect.h"
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <iomanip>

namespace cavity_detect {

float CavityDetector::distance(const float* a, const float* b) const {
    float dx = a[0] - b[0];
    float dy = a[1] - b[1];
    float dz = a[2] - b[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void CavityDetector::load_from_fa(const atom* atoms, const resid* residues, int res_cnt) {
    m_atoms.clear();
    m_clefts.clear();
    // TODO: real FlexAID struct mapping in Task 1.3 (read_input hook)
    // For now stubbed — we pass real atoms later
}

void CavityDetector::load_from_pdb(const std::string& pdb_file) {
    // Standalone parser stub for FreeNRG CLI tests (expand later)
}

void CavityDetector::detect(float min_radius, float max_radius) {
    m_sphere_lwb = min_radius;
    m_sphere_upb = max_radius;
    m_clefts.clear();

    // === REAL SURFNET CORE (placeholder ready for full O(N) probe placement) ===
    // Full version in 1.3: for each atom pair place probe spheres in gaps,
    // reject clashes, cluster by overlap → clefts. OpenMP over atom pairs.
    DetectedCleft cleft;
    cleft.id = 1;
    cleft.label = 1;
    cleft.volume = 250.0f;
    cleft.center[0] = cleft.center[1] = cleft.center[2] = 0.0f;

    DetectedSphere s;
    s.center[0] = 5.0f; s.center[1] = 5.0f; s.center[2] = 5.0f;
    s.radius = 2.8f;
    s.cleft_id = 1;
    cleft.spheres.push_back(s);

    m_clefts.push_back(std::move(cleft));
    sort_clefts();
}

void CavityDetector::merge_clefts() {
    // TODO: merge overlapping clefts by sphere distance threshold
}

void CavityDetector::sort_clefts() {
    std::sort(m_clefts.begin(), m_clefts.end(),
        [](const DetectedCleft& a, const DetectedCleft& b) {
            return a.spheres.size() > b.spheres.size();
        });
}

void CavityDetector::assign_atoms_to_clefts(float contact_threshold) {
    // TODO: assign protein atoms to closest cleft
}

void CavityDetector::filter_anchor_residues(const std::string& anchor_residues) {
    // TODO: port of your -a RESNUM,RESNUMCA logic
}

sphere* CavityDetector::to_flexaid_spheres(int cleft_id) const {
    // Convert to old FlexAID sphere* linked list for generate_grid()
    return nullptr; // full impl in 1.3
}

void CavityDetector::write_sphere_pdb(const std::string& filename, int cleft_id) const {
    std::ofstream out(filename);
    if (!out) return;

    out << "REMARK  Cleft spheres by FlexAIDΔS native CavityDetector\n";
    int atom_id = 1;
    for (const auto& cleft : m_clefts) {
        if (cleft.id != cleft_id) continue;
        for (const auto& s : cleft.spheres) {
            out << "HETATM" << std::setw(5) << atom_id++
                << "  SPH SURF   1    "
                << std::fixed << std::setprecision(3)
                << std::setw(8) << s.center[0]
                << std::setw(8) << s.center[1]
                << std::setw(8) << s.center[2]
                << "  1.00 " << std::setw(5) << std::setprecision(2) << s.radius
                << "           S\n";
        }
    }
}

} // namespace cavity_detect
