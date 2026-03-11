// statmech.cpp — Statistical Mechanics Engine implementation
//
// Notation:
//   β  = 1/(kB T)
//   Z  = Σ_i  n_i exp(−β E_i)        (canonical partition function)
//   F  = −kT ln Z                      (Helmholtz free energy)
//   ⟨E⟩ = (1/Z) Σ_i  n_i E_i exp(−β E_i)
//   C_v = (⟨E²⟩ − ⟨E⟩²) / (kT²)      (heat capacity)
//   S  = (⟨E⟩ − F) / T                 (entropy)
//
// All sums use log-sum-exp for numerical stability when energies span
// hundreds of kcal/mol (common in docking).
//
// Hardware acceleration dispatch (compile-time):
//   1. AVX-512  (__AVX512F__)  — 8 doubles/cycle vectorised exp/reduction
//   2. Eigen    (FLEXAIDS_HAS_EIGEN) — vectorised array ops
//   3. OpenMP   (_OPENMP)     — parallel histogram binning in WHAM
//   4. Scalar   (always available)

#include "statmech.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <stdexcept>

#ifdef __AVX512F__
#  include <immintrin.h>
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

#ifdef FLEXAIDS_HAS_EIGEN
#  include <Eigen/Dense>
#endif

namespace statmech {

// ─── construction ────────────────────────────────────────────────────────────

StatMechEngine::StatMechEngine(double temperature_K)
    : T_(temperature_K)
    , beta_(1.0 / (kB_kcal * temperature_K))
{
    if (temperature_K <= 0.0)
        throw std::invalid_argument("StatMechEngine: temperature must be > 0");
}

// ─── add_sample ──────────────────────────────────────────────────────────────

void StatMechEngine::add_sample(double energy, int multiplicity) {
    ensemble_.push_back({energy, multiplicity});
}

// ─── log_sum_exp ─────────────────────────────────────────────────────────────
// Numerically stable log(Σ exp(x_i)) with multi-backend acceleration.

double StatMechEngine::log_sum_exp(std::span<const double> x) {
    if (x.empty()) return -1e308;

    double x_max = *std::max_element(x.begin(), x.end());
    if (x_max <= -1e308) return x_max;

    const std::size_t N = x.size();
    double sum = 0.0;

#ifdef __AVX512F__
    // AVX-512: process 8 doubles per cycle with vectorised exp
    {
        __m512d vmax = _mm512_set1_pd(x_max);
        __m512d vacc = _mm512_setzero_pd();
        std::size_t i = 0;
        for (; i + 7 < N; i += 8) {
            __m512d vx   = _mm512_loadu_pd(x.data() + i);
            __m512d vdif = _mm512_sub_pd(vx, vmax);
            // _mm512_exp_pd requires SVML; available with Intel compilers & GCC+SVML
            __m512d vexp = _mm512_exp_pd(vdif);
            vacc = _mm512_add_pd(vacc, vexp);
        }
        sum = _mm512_reduce_add_pd(vacc);
        // Scalar tail
        for (; i < N; ++i)
            sum += std::exp(x[i] - x_max);
    }
#elif defined(FLEXAIDS_HAS_EIGEN)
    // Eigen: vectorised exp over the shifted array
    {
        Eigen::Map<const Eigen::ArrayXd> xarr(x.data(), static_cast<Eigen::Index>(N));
        sum = (xarr - x_max).exp().sum();
    }
#else
    // Scalar fallback
    for (std::size_t i = 0; i < N; ++i)
        sum += std::exp(x[i] - x_max);
#endif

    return x_max + std::log(sum);
}

// ─── compute ─────────────────────────────────────────────────────────────────

Thermodynamics StatMechEngine::compute() const {
    if (ensemble_.empty())
        throw std::runtime_error("StatMechEngine::compute: empty ensemble");

    const std::size_t N = ensemble_.size();

    // Build array of log-weights:  w_i = ln(n_i) − β E_i
    std::vector<double> log_w(N);

#ifdef FLEXAIDS_HAS_EIGEN
    // Eigen-vectorised log-weight construction
    {
        Eigen::ArrayXd counts(static_cast<Eigen::Index>(N));
        Eigen::ArrayXd energies(static_cast<Eigen::Index>(N));
        for (std::size_t i = 0; i < N; ++i) {
            counts(static_cast<Eigen::Index>(i))   = static_cast<double>(ensemble_[i].count);
            energies(static_cast<Eigen::Index>(i)) = ensemble_[i].energy;
        }
        Eigen::ArrayXd lw = counts.log() - beta_ * energies;
        Eigen::Map<Eigen::ArrayXd>(log_w.data(), static_cast<Eigen::Index>(N)) = lw;
    }
#else
    for (std::size_t i = 0; i < N; ++i)
        log_w[i] = std::log(static_cast<double>(ensemble_[i].count)) -
                   beta_ * ensemble_[i].energy;
#endif

    double lnZ = log_sum_exp(log_w);

    // ⟨E⟩ and ⟨E²⟩ via Boltzmann probabilities
    double E_avg  = 0.0;
    double E2_avg = 0.0;

#ifdef FLEXAIDS_HAS_EIGEN
    {
        Eigen::Map<const Eigen::ArrayXd> lw_arr(log_w.data(), static_cast<Eigen::Index>(N));
        Eigen::ArrayXd p = (lw_arr - lnZ).exp();  // Boltzmann probabilities

        Eigen::ArrayXd E(static_cast<Eigen::Index>(N));
        for (std::size_t i = 0; i < N; ++i)
            E(static_cast<Eigen::Index>(i)) = ensemble_[i].energy;

        E_avg  = (p * E).sum();
        E2_avg = (p * E * E).sum();
    }
#else
    for (std::size_t i = 0; i < N; ++i) {
        double p_i = std::exp(log_w[i] - lnZ);
        double Ei  = ensemble_[i].energy;
        E_avg  += p_i * Ei;
        E2_avg += p_i * Ei * Ei;
    }
#endif

    double kT  = kB_kcal * T_;
    double var = E2_avg - E_avg * E_avg;

    Thermodynamics th;
    th.temperature    = T_;
    th.log_Z          = lnZ;
    th.free_energy    = -kT * lnZ;
    th.mean_energy    = E_avg;
    th.mean_energy_sq = E2_avg;
    th.heat_capacity  = var / (kT * kT);
    th.entropy        = (E_avg - th.free_energy) / T_;
    th.std_energy     = std::sqrt(std::max(0.0, var));
    return th;
}

// ─── boltzmann_weights ───────────────────────────────────────────────────────

std::vector<double> StatMechEngine::boltzmann_weights() const {
    if (ensemble_.empty()) return {};

    const std::size_t N = ensemble_.size();
    std::vector<double> log_w(N);

#ifdef FLEXAIDS_HAS_EIGEN
    {
        Eigen::ArrayXd counts(static_cast<Eigen::Index>(N));
        Eigen::ArrayXd energies(static_cast<Eigen::Index>(N));
        for (std::size_t i = 0; i < N; ++i) {
            counts(static_cast<Eigen::Index>(i))   = static_cast<double>(ensemble_[i].count);
            energies(static_cast<Eigen::Index>(i)) = ensemble_[i].energy;
        }
        Eigen::ArrayXd lw = counts.log() - beta_ * energies;
        Eigen::Map<Eigen::ArrayXd>(log_w.data(), static_cast<Eigen::Index>(N)) = lw;
    }
#else
    for (std::size_t i = 0; i < N; ++i)
        log_w[i] = std::log(static_cast<double>(ensemble_[i].count)) -
                   beta_ * ensemble_[i].energy;
#endif

    double lnZ = log_sum_exp(log_w);

    std::vector<double> w(N);

#ifdef FLEXAIDS_HAS_EIGEN
    {
        Eigen::Map<const Eigen::ArrayXd> lw_arr(log_w.data(), static_cast<Eigen::Index>(N));
        Eigen::Map<Eigen::ArrayXd>(w.data(), static_cast<Eigen::Index>(N)) =
            (lw_arr - lnZ).exp();
    }
#else
    for (std::size_t i = 0; i < N; ++i)
        w[i] = std::exp(log_w[i] - lnZ);
#endif

    return w;
}

// ─── delta_G ─────────────────────────────────────────────────────────────────
// ΔG = F_this − F_ref = −kT (ln Z_this − ln Z_ref)

double StatMechEngine::delta_G(const StatMechEngine& reference) const {
    auto this_th = this->compute();
    auto ref_th  = reference.compute();
    double kT = kB_kcal * T_;
    return -kT * (this_th.log_Z - ref_th.log_Z);
}

// ─── Helmholtz convenience ───────────────────────────────────────────────────

double StatMechEngine::helmholtz(std::span<const double> energies, double T) {
    if (energies.empty())
        throw std::invalid_argument("helmholtz: empty energy list");

    const std::size_t N = energies.size();
    double beta = 1.0 / (kB_kcal * T);
    std::vector<double> neg_beta_E(N);

#ifdef FLEXAIDS_HAS_EIGEN
    {
        Eigen::Map<const Eigen::ArrayXd> E(energies.data(), static_cast<Eigen::Index>(N));
        Eigen::Map<Eigen::ArrayXd>(neg_beta_E.data(), static_cast<Eigen::Index>(N)) =
            -beta * E;
    }
#else
    for (std::size_t i = 0; i < N; ++i)
        neg_beta_E[i] = -beta * energies[i];
#endif

    double lnZ = log_sum_exp(neg_beta_E);
    return -(kB_kcal * T) * lnZ;
}

// ─── replica exchange ────────────────────────────────────────────────────────

std::vector<Replica>
StatMechEngine::init_replicas(std::span<const double> temperatures) {
    std::vector<Replica> reps;
    reps.reserve(temperatures.size());
    int id = 0;
    for (double T : temperatures) {
        Replica r;
        r.id             = id++;
        r.temperature    = T;
        r.beta           = 1.0 / (kB_kcal * T);
        r.current_energy = 0.0;
        reps.push_back(r);
    }
    return reps;
}

bool StatMechEngine::attempt_swap(Replica& a, Replica& b, std::mt19937& rng) {
    // Metropolis criterion:
    //   Δ = (β_a − β_b)(E_a − E_b)
    //   P_accept = min(1, exp(Δ))
    double delta = (a.beta - b.beta) * (a.current_energy - b.current_energy);
    if (delta >= 0.0) {
        std::swap(a.current_energy, b.current_energy);
        return true;
    }
    std::uniform_real_distribution<double> U(0.0, 1.0);
    if (U(rng) < std::exp(delta)) {
        std::swap(a.current_energy, b.current_energy);
        return true;
    }
    return false;
}

// ─── WHAM ────────────────────────────────────────────────────────────────────
// Weighted Histogram Analysis Method (Kumar et al. 1992)
// Simplified single-window version for post-hoc reweighting of GA ensemble.
//
// Acceleration:
//   - AVX-512 + OpenMP: vectorised bin index computation + parallel private histograms
//   - OpenMP only: parallel for with private histogram reduction
//   - Scalar: sequential fallback

std::vector<WHAMBin> StatMechEngine::wham(
    std::span<const double> energies,
    std::span<const double> coordinates,
    double temperature,
    int    n_bins,
    int    max_iter,
    double tolerance)
{
    if (energies.size() != coordinates.size())
        throw std::invalid_argument("wham: energies and coordinates size mismatch");
    if (energies.empty() || n_bins <= 0)
        throw std::invalid_argument("wham: invalid input");

    const std::size_t N = energies.size();
    const int Ni = static_cast<int>(N);
    double beta = 1.0 / (kB_kcal * temperature);

    // Find coordinate range
    double cmin = *std::min_element(coordinates.begin(), coordinates.end());
    double cmax = *std::max_element(coordinates.begin(), coordinates.end());
    double bin_w = (cmax - cmin) / n_bins;
    if (bin_w <= 0.0) bin_w = 1.0;
    double inv_bw = 1.0 / bin_w;

    // Histogram + Boltzmann-weighted histogram
    std::vector<double> raw_count(static_cast<std::size_t>(n_bins), 0.0);
    std::vector<double> boltz_sum(static_cast<std::size_t>(n_bins), 0.0);

    // ── Histogram binning with hardware dispatch ──

#if defined(__AVX512F__) && defined(_OPENMP)
    // AVX-512 bin-index computation + OpenMP parallel private histograms
    {
        int n_threads = omp_get_max_threads();
        std::vector<std::vector<double>> t_raw(n_threads,
            std::vector<double>(static_cast<std::size_t>(n_bins), 0.0));
        std::vector<std::vector<double>> t_boltz(n_threads,
            std::vector<double>(static_cast<std::size_t>(n_bins), 0.0));

        __m512d vcmin  = _mm512_set1_pd(cmin);
        __m512d vinvbw = _mm512_set1_pd(inv_bw);
        __m512d vnbeta = _mm512_set1_pd(-beta);

        #pragma omp parallel
        {
            int tid   = omp_get_thread_num();
            int chunk = (Ni + n_threads - 1) / n_threads;
            int start = tid * chunk;
            int end   = std::min(start + chunk, Ni);

            int i = start;
            // AVX-512 vectorised bin index computation (8 doubles at a time)
            for (; i + 7 < end; i += 8) {
                __m512d vc = _mm512_loadu_pd(coordinates.data() + i);
                __m512d ve = _mm512_loadu_pd(energies.data() + i);

                // bin = (coord - cmin) * inv_bw
                __m512d vrel = _mm512_mul_pd(_mm512_sub_pd(vc, vcmin), vinvbw);
                __m256i vbin = _mm512_cvttpd_epi32(vrel);

                // Boltzmann factor: exp(-beta * E)
                __m512d vbf = _mm512_exp_pd(_mm512_mul_pd(vnbeta, ve));

                alignas(32) int bins[8];
                alignas(64) double bfs[8];
                _mm256_storeu_si256((__m256i*)bins, vbin);
                _mm512_storeu_pd(bfs, vbf);

                for (int k = 0; k < 8; ++k) {
                    int b = std::min(std::max(bins[k], 0), n_bins - 1);
                    t_raw[tid][static_cast<std::size_t>(b)]  += 1.0;
                    t_boltz[tid][static_cast<std::size_t>(b)] += bfs[k];
                }
            }
            // Scalar tail
            for (; i < end; ++i) {
                int b = static_cast<int>((coordinates[i] - cmin) * inv_bw);
                b = std::min(std::max(b, 0), n_bins - 1);
                t_raw[tid][static_cast<std::size_t>(b)]  += 1.0;
                t_boltz[tid][static_cast<std::size_t>(b)] += std::exp(-beta * energies[i]);
            }
        }
        // Reduce thread-private histograms
        for (auto& tr : t_raw)
            for (int b = 0; b < n_bins; ++b)
                raw_count[static_cast<std::size_t>(b)] += tr[static_cast<std::size_t>(b)];
        for (auto& tb : t_boltz)
            for (int b = 0; b < n_bins; ++b)
                boltz_sum[static_cast<std::size_t>(b)] += tb[static_cast<std::size_t>(b)];
    }

#elif defined(_OPENMP)
    // OpenMP parallel histogram with per-thread private bins
    {
        int n_threads = omp_get_max_threads();
        std::vector<std::vector<double>> t_raw(n_threads,
            std::vector<double>(static_cast<std::size_t>(n_bins), 0.0));
        std::vector<std::vector<double>> t_boltz(n_threads,
            std::vector<double>(static_cast<std::size_t>(n_bins), 0.0));

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < Ni; ++i) {
            int tid = omp_get_thread_num();
            int b = static_cast<int>((coordinates[i] - cmin) * inv_bw);
            b = std::min(std::max(b, 0), n_bins - 1);
            t_raw[tid][static_cast<std::size_t>(b)]  += 1.0;
            t_boltz[tid][static_cast<std::size_t>(b)] += std::exp(-beta * energies[i]);
        }
        for (auto& tr : t_raw)
            for (int b = 0; b < n_bins; ++b)
                raw_count[static_cast<std::size_t>(b)] += tr[static_cast<std::size_t>(b)];
        for (auto& tb : t_boltz)
            for (int b = 0; b < n_bins; ++b)
                boltz_sum[static_cast<std::size_t>(b)] += tb[static_cast<std::size_t>(b)];
    }

#else
    // Scalar fallback
    for (std::size_t i = 0; i < N; ++i) {
        int b = static_cast<int>((coordinates[i] - cmin) / bin_w);
        if (b < 0) b = 0;
        if (b >= n_bins) b = n_bins - 1;
        raw_count[static_cast<std::size_t>(b)] += 1.0;
        boltz_sum[static_cast<std::size_t>(b)] += std::exp(-beta * energies[i]);
    }
#endif

    // ── Free energy self-consistency iteration ──
    // Eigen-accelerated when available for vectorised log/subtract/abs

    std::vector<double> f_old(static_cast<std::size_t>(n_bins), 0.0);
    std::vector<double> f_new(static_cast<std::size_t>(n_bins), 0.0);

    for (int iter = 0; iter < max_iter; ++iter) {

#ifdef FLEXAIDS_HAS_EIGEN
        {
            Eigen::Map<const Eigen::ArrayXd> rc(raw_count.data(), n_bins);
            Eigen::Map<const Eigen::ArrayXd> bs(boltz_sum.data(), n_bins);
            Eigen::Map<Eigen::ArrayXd> fn(f_new.data(), n_bins);

            // F_b = -kT * ln(boltz_sum / raw_count)  where raw_count > 0
            Eigen::ArrayXd mask = (rc > 0.0).cast<double>();
            Eigen::ArrayXd safe_rc = (rc > 0.0).select(rc, Eigen::ArrayXd::Ones(n_bins));
            fn = mask * (-(kB_kcal * temperature) * (bs / safe_rc).log());
        }
#else
        for (int b = 0; b < n_bins; ++b) {
            if (raw_count[static_cast<std::size_t>(b)] > 0.0) {
                f_new[static_cast<std::size_t>(b)] = -(kB_kcal * temperature) *
                    std::log(boltz_sum[static_cast<std::size_t>(b)] /
                             raw_count[static_cast<std::size_t>(b)]);
            } else {
                f_new[static_cast<std::size_t>(b)] = 0.0;
            }
        }
#endif

        // Shift so minimum = 0
        double fmin = *std::min_element(f_new.begin(), f_new.end());
        for (auto& f : f_new) f -= fmin;

        // Check convergence
        double maxdiff = 0.0;
#ifdef FLEXAIDS_HAS_EIGEN
        {
            Eigen::Map<const Eigen::ArrayXd> fn(f_new.data(), n_bins);
            Eigen::Map<const Eigen::ArrayXd> fo(f_old.data(), n_bins);
            maxdiff = (fn - fo).abs().maxCoeff();
        }
#else
        for (int b = 0; b < n_bins; ++b)
            maxdiff = std::max(maxdiff,
                std::abs(f_new[static_cast<std::size_t>(b)] -
                         f_old[static_cast<std::size_t>(b)]));
#endif
        f_old = f_new;
        if (maxdiff < tolerance) break;
    }

    // Build output
    std::vector<WHAMBin> result(static_cast<std::size_t>(n_bins));
    for (int b = 0; b < n_bins; ++b) {
        result[static_cast<std::size_t>(b)].coord_center = cmin + (b + 0.5) * bin_w;
        result[static_cast<std::size_t>(b)].count        = raw_count[static_cast<std::size_t>(b)];
        result[static_cast<std::size_t>(b)].free_energy  = f_new[static_cast<std::size_t>(b)];
    }
    return result;
}

// ─── thermodynamic integration ───────────────────────────────────────────────
// ΔG = ∫₀¹ ⟨∂V/∂λ⟩_λ dλ   (trapezoidal rule)

double StatMechEngine::thermodynamic_integration(std::span<const TIPoint> points) {
    if (points.size() < 2)
        throw std::invalid_argument("TI requires at least 2 points");

    double integral = 0.0;
    for (std::size_t i = 1; i < points.size(); ++i) {
        double dl = points[i].lambda - points[i-1].lambda;
        integral += 0.5 * dl * (points[i].dV_dlambda + points[i-1].dV_dlambda);
    }
    return integral;
}

// ─── BoltzmannLUT ────────────────────────────────────────────────────────────

BoltzmannLUT::BoltzmannLUT(double beta, double e_min, double e_max, int n_bins)
    : beta_(beta)
    , e_min_(e_min)
    , n_bins_(n_bins)
    , table_(static_cast<std::size_t>(n_bins))
{
    double range = e_max - e_min;
    if (range <= 0.0) range = 1.0;
    inv_bin_width_ = n_bins / range;

#ifdef FLEXAIDS_HAS_EIGEN
    // Eigen-vectorised LUT initialisation
    {
        Eigen::ArrayXd idx = Eigen::ArrayXd::LinSpaced(n_bins, 0.5, n_bins - 0.5);
        Eigen::ArrayXd E = e_min + idx * (range / n_bins);
        Eigen::Map<Eigen::ArrayXd>(table_.data(), n_bins) = (-beta * E).exp();
    }
#else
    for (int i = 0; i < n_bins; ++i) {
        double e = e_min + (static_cast<double>(i) + 0.5) * range / n_bins;
        table_[static_cast<std::size_t>(i)] = std::exp(-beta * e);
    }
#endif
}

double BoltzmannLUT::operator()(double energy) const noexcept {
    int idx = static_cast<int>((energy - e_min_) * inv_bin_width_);
    if (idx < 0) idx = 0;
    if (idx >= n_bins_) idx = n_bins_ - 1;
    return table_[static_cast<std::size_t>(idx)];
}

}  // namespace statmech
