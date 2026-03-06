// shannon_cuda.cuh — CUDA Shannon histogram kernel declarations
#pragma once

#ifdef FLEXAIDS_USE_CUDA
#include <cuda_runtime.h>

// Opaque GPU context for Shannon histogram
struct ShannonCudaCtx {
    double* d_energies = nullptr;
    int*    d_bins     = nullptr;
    int     capacity   = 0;   // max energies allocated
    int     num_bins   = 0;
};

// Allocate device buffers for up to `max_n` energy values and `num_bins` bins.
void shannon_cuda_init(ShannonCudaCtx& ctx, int max_n, int num_bins);

// Free device buffers.
void shannon_cuda_shutdown(ShannonCudaCtx& ctx);

// Copy `n` energies to device, run histogram kernel, copy bin counts back.
// `bins_out` must be caller-allocated with ctx.num_bins ints.
void shannon_cuda_histogram(ShannonCudaCtx& ctx,
                             const double*   energies_host,
                             int             n,
                             double          min_v,
                             double          bin_width,
                             int*            bins_out);

#endif // FLEXAIDS_USE_CUDA
