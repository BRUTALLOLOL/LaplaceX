#pragma once

#include <vector>
#include <complex>
#include <string>

namespace laplacex {

// ---------------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------------

/**
 * A single exponential component recovered from a signal.
 *
 *   contribution(t) = amplitude * exp(-rate * t)
 *
 * For oscillatory modes, frequency != 0 and the contribution is
 *   amplitude * exp(-rate * t) * cos(frequency * t + phase)
 */
struct ExponentialComponent {
    double amplitude   = 0.0;  ///< A_k
    double rate        = 0.0;  ///< lambda_k  (decay rate, >= 0)
    double decay_time  = 0.0;  ///< 1 / lambda_k  (relaxation time)
    double frequency   = 0.0;  ///< omega_k (0 for purely real modes)
    double phase       = 0.0;  ///< phi_k
    double significance = 0.0; ///< normalized absolute amplitude (0–1)
    double log_time    = 0.0;  ///< x_k = ln(t_k) – position in log plane
};

// ---------------------------------------------------------------------------
// Critical-points decomposition
// ---------------------------------------------------------------------------

struct DecomposeConfig {
    // Wavelet
    std::string wavelet    = "mexh";  ///< only Mexican Hat is implemented
    int    n_scales        = 16;      ///< number of logarithmically spaced scales
    double scale_min       = 0.05;    ///< smallest wavelet scale (log-space units)
    double scale_max       = 0.5;     ///< largest  wavelet scale (log-space units)

    // Log-space resampling
    int    n_log_points    = 256;     ///< uniform grid size in x = ln(t)

    // Skeleton / ridge tracking
    double ridge_threshold = 0.05;   ///< min |W| / max|W| to keep a ridge line
    int    min_ridge_len   = 4;      ///< minimum ridge length (in scale steps)

    // Noise & filtering
    double noise_sigma     = -1.0;   ///< <0 → MAD estimate from local interpolation residuals
    double snr_threshold   = 3.0;    ///< min amplitude / noise to keep component

    // Model selection
    bool   use_aic         = true;   ///< legacy enable switch for information-criterion pruning
    std::string model_selection = "bic"; ///< "bic" (parsimonious) or "aicc"

    // Advanced
    /// Legacy name: skip joint nonlinear rate/amplitude refinement.
    /// Linear amplitudes are still fitted by QR; raw ridge rates are approximate.
    bool   disable_als     = false;
};

/**
 * Fit a real sum of decaying exponentials, with signed amplitudes allowed.
 * Mexican-Hat zero crossings initialize a joint QR/LM fit in log time.
 * This routine does not estimate oscillation frequencies.
 *
 * @param t       Strictly increasing time values (t > 0)
 * @param f       Signal values at t
 * @param cfg     Algorithm configuration
 * @return        Sorted (by significance, descending) list of components
 */
std::vector<ExponentialComponent> decompose(
    const std::vector<double>& t,
    const std::vector<double>& f,
    const DecomposeConfig& cfg = {},
    int max_components = 5
);

/**
 * Fit using log-spaced rate-grid proposals instead of wavelet proposals.
 * Shares the same joint QR/LM optimizer, noise tests and model selection.
 *
 * @param t              Strictly increasing time values (t > 0)
 * @param f              Signal values at t
 * @param cfg            Decomposition configuration
 * @param max_components Maximum number of components to extract (default 5)
 * @param energy_tol     Relative residual stopping target when model selection
 *                       is disabled (default 1e-6)
 * @return               Sorted (by significance, descending) list of components
 */
std::vector<ExponentialComponent> decompose_deflate(
    const std::vector<double>& t,
    const std::vector<double>& f,
    const DecomposeConfig& cfg = {},
    int max_components = 5,
    double energy_tol = 1e-6
);

// ---------------------------------------------------------------------------
// Reconstruction
// ---------------------------------------------------------------------------

/**
 * Reconstruct f(t) from a list of exponential components.
 */
std::vector<double> reconstruct(
    const std::vector<ExponentialComponent>& components,
    const std::vector<double>& t_eval
);

/**
 * Compute the Laplace image F(s) analytically from components:
 *   F(s) = sum_k  A_k / (s + lambda_k)          (continuous exponentials)
 */
std::vector<std::complex<double>> reconstruct_laplace(
    const std::vector<ExponentialComponent>& components,
    const std::vector<std::complex<double>>& s_eval
);

} // namespace laplacex
