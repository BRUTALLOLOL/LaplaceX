# Log-plane wavelet recovery and Python comparisons

Measured on 20 October 2026 after the log-wavelet rewrite. Times are medians
of five runs after one warmup, with shuffled method order and one BLAS thread.

![Wavelet recovery comparison](comparison/wavelet_comparison.png)

## Component recovery

| Signal | LaplaceX default, ms | Fastest passing SciPy auto-K, ms | Speed ratio | Modes: LaplaceX / truth |
|---|---:|---:|---:|---:|
| single_uniform | 0.49 | 3.03 | 6.14× | 1 / 1 |
| two_uniform | 4.03 | 12.45 | 3.09× | 2 / 2 |
| three_log_pilot | 2.87 | 37.75 | 13.14× | 3 / 3 |
| two_log_pilot | 2.41 | 7.31 | 3.03× | 2 / 2 |
| three_log_holdout | 2.86 | 17.69 | 6.18× | 3 / 3 |
| close_rates_holdout | 12.87 | 54.22 | 4.21× | 2 / 2 |
| weak_fast | 4.25 | 18.88 | 4.45× | 3 / 3 |
| weak_middle | 3.76 | 13.57 | 3.61× | 3 / 3 |
| weak_slow | 4.32 | 42.02 | 9.73× | 3 / 3 |
| noisy_log_seed7 | 4.76 | 39.32 | 8.26× | 3 / 3 |
| noisy_log_seed19 | 4.38 | 37.60 | 8.58× | 3 / 3 |
| noisy_log_seed41 | 4.85 | 36.93 | 7.62× | 3 / 3 |

A passing automatic-count fit must recover the true count. Clean cases also
require normalized max reconstruction error ≤1e-8 and rate error ≤1e-5.
Noisy cases require clean-signal RMSE ≤noise σ and rate error ≤10%.
The comparator is the fastest passing automatic-K SciPy configuration among
the original residual-stop TRF fitter and TRF/LM BIC fitters. All outputs,
including failures and incorrect counts, remain in the full tables below.

## What changed

- Analyze `f(exp(x))`, rather than `exp(x)*f(exp(x))`. The negative log
  derivative of one exponential has the profile `A*exp(z-exp(z))`,
  `z=x+log(rate)`: rate changes translate the profile without changing
  its height. The old time factor favored slow components.
- Track persistent Mexican-Hat zero crossings from fine to coarse scales,
  retaining branches when coarse scales merge them. Initialize all persistent
  components in one fit instead of repeatedly fitting incomplete models.
- Fit rates and signed amplitudes jointly with an analytic Jacobian, damped
  Levenberg–Marquardt steps and scaled Householder QR. Final fitting uses
  the original observations, so log-grid interpolation is initialization only.
- Propagate sample noise through interpolation and wavelet filter weights,
  including covariance from reused observations. Condition residual proposals
  on both the existing amplitudes and rates.
- Use BIC by default, refit candidate deletions, and screen nearly dependent
  modes using conditional amplitude uncertainty. AICc remains configurable.
- Respect the component budget in both stages. Handle constants, noise-only
  signals, and large/small amplitude scales without spurious decay modes or
  overflowing the Python residual-energy calculation.

- Share an exact single-exponential shortcut between both routes. It accepts
  a model only after fitting every sample at the common numerical floor.
  Its single-mode speed benefit is not attributed to wavelet analysis.

## Accuracy of automatic model selection

| Weak mode | LaplaceX modes / max abs error | SciPy residual-stop modes / max abs error | SciPy BIC / known-K accuracy |
|---|---:|---:|---|
| weak_fast | 3 / 4.44e-16 | 2 / 4.25e-04 | Also recovers 3 modes; error down to 6.66e-16 |
| weak_middle | 3 / 4.44e-16 | 2 / 3.79e-04 | Also recovers 3 modes; error down to 4.44e-16 |
| weak_slow | 3 / 5.20e-16 | 2 / 8.53e-04 | Also recovers 3 modes; error down to 4.44e-16 |

The 1000:1 amplitude cases expose the old SciPy automatic baseline's
`1e-6` residual-energy stopping rule: it can accept two modes before a weak
third mode is recovered. LaplaceX finds all three. This is an advantage over
that automatic stopping policy, not universal numerical accuracy superiority:
the stronger SciPy BIC and known-K variants also recover the modes accurately.
Float64 roundoff differences between successful fits are not meaningful accuracy wins.

## Independent validation

The separately seeded validation recovered the correct count in **160/160** mixtures
with the rate-grid fallback disabled. It includes 80 clean and 80 noisy signals,
1–4 modes, linear and logarithmic sampling, and independently varied amplitudes/rates.
Maximum clean rate error: `1.784e-12` relative.
For noise σ=0.001, median / 95th-percentile / maximum rate errors were
**0.26% / 4.55% / 7.30%**.
Correct counts do not imply exact rates in noisy, poorly conditioned mixtures.
These are finite synthetic validation results, not a guarantee for all signals.

[Full validation protocol and individual outcomes](comparison/validation.json).
Model-selection development used seed 20261023; the separate validation
uses seed 20261029. Later sparse-variance and exact-single-mode performance
optimizations were revalidated on the same protocol; no selection parameters
were tuned to the second seed.
Validation timings are single observations and are not used for speed claims.

## Wavelet ablation

`LaplaceX no wavelets` bypasses CWT and searches a 64-point log-rate grid.
It shares the new QR/LM solver, noise estimation, conditional candidate test,
BIC and deletion refits. This isolates initialization and candidate generation
from the optimizer rewrite. It removes the entire log-CWT stage; it does not
separately compare linear-time and log-time wavelets.

The default public route also includes Python fallback/noise checks. The
native no-wavelet entry omits that routing overhead; the explicitly
wavelet-only rows omit fallback checks and provide an additional comparison.

With the same optimizer, default wavelet initialization is **8.91×** faster
on the widely separated three-mode case and **3.62×** faster on its holdout.
It is slower on the two-uniform and close-rate cases (ratios **0.90×** and
**0.67×**). These ratios compare the native rate-grid route
with the public wavelet route, including its Python routing overhead.

![Inverse speed and wavelet ablation](comparison/comparison.png)

The default configuration now uses 16 scales over [0.05,0.5] and 256 log points.
The older fine preset remains in the matrix. Both automatic fallback behavior
and the explicitly wavelet-only route are tested. Excess residual consistent
with estimated noise no longer triggers a redundant rate-grid fit.

## Before and after on identical fixtures

These compare the previous default library with the rewritten default library;
they include the combined effects of wavelet recovery, fitting and selection.
The ablation above is the narrower test of wavelet contribution.

| Case | Previous default ms / modes | New default ms / modes | Max abs error before → after |
|---|---:|---:|---:|
| single_uniform | 25.32 / 1 | 0.49 / 1 | 7.26e-14 → 3.33e-16 |
| two_uniform | 294.21 / 2 | 4.03 / 2 | 3.36e-12 → 4.44e-16 |
| three_log_pilot | 628.57 / 3 | 2.87 / 3 | 2.54e-12 → 9.99e-16 |
| two_log_pilot | 113.84 / 2 | 2.41 / 2 | 1.19e-13 → 2.22e-16 |
| three_log_holdout | 646.35 / 3 | 2.86 / 3 | 5.26e-12 → 4.44e-16 |
| close_rates_holdout | 1004.65 / 2 | 12.87 / 2 | 8.62e-05 → 4.44e-16 |
| noisy_log_seed7 | 574.99 / 3 | 4.76 / 3 | 3.74e-04 → 3.74e-04 |
| noisy_log_seed19 | 581.83 / 3 | 4.38 / 3 | 1.58e-04 → 1.58e-04 |
| noisy_log_seed41 | 582.79 / 3 | 4.85 / 3 | 1.77e-04 → 1.77e-04 |

Historical raw files and benchmark scripts are preserved in
[comparison-before-wavelet](comparison-before-wavelet/decomposition.json).

## Forward and inverse comparisons

| Inverse workload (64 times) | LaplaceX, ms / max error | Fastest passing mpmath, ms / max error | Speed ratio |
|---|---:|---:|---:|
| exp | 1.54 / 1.61e-14 (LaplaceX N64) | 159.59 / 0.00e+00 (mpmath Cohen dps30 auto) | 103.6× |
| repeated_pole | 1.71 / 3.16e-15 (LaplaceX N64) | 112.05 / 3.76e-10 (mpmath Cohen degree32) | 65.5× |
| damped_cos | 4.78 / 3.32e-12 (LaplaceX N64) | 270.64 / 0.00e+00 (mpmath Cohen dps30 auto) | 56.7× |
| sin10 | 91.56 / 1.50e-14 (LaplaceX N512) | 194.42 / 1.11e-16 (mpmath Cohen dps30 auto) | 2.1× |
| diffusion | No tested setting passed | 134.37 / 1.29e-10 (mpmath Cohen degree32) | — |
| delayed_step | No tested setting passed | 15644.64 / 1.07e-13 (mpmath deHoog degree64) | — |

Inverse ratios select the fastest passing tested setting in each package.
Every point must satisfy `abs(error) <= 1e-10+1e-8*abs(reference)`.
The Weeks inverse and forward quadrature do not use wavelet decomposition.
No wavelet-related accuracy gain is attributed to those paths. The inverse
still has limitations on singular/discontinuous signals. SciPy's rational
impulse response has richer inputs and is retained as a separate comparison.

## Fairness, settings, and limits

- Python `3.14.7 (tags/v3.14.7:823f032, Aug  5 2026, 10:51:32) [MSC v.1944 64 bit (AMD64)]`; `Windows-10-10.0.19045-SP0`.
- NumPy `2.5.3`, SciPy `1.18.1`, mpmath `1.3.0`
  with its `python` backend (gmpy2 not installed).
- LaplaceX: Release Windows x64 wheel, Zig 0.16.0 / Clang 21.1.0, C++17, -O3.
  Extension and script SHA-256 hashes, all run times and warnings are in JSON.
- Sections run serially. Setup, optimization and output conversion are timed;
  imports and reference generation are not. Medians are local observations,
  not cross-machine performance guarantees or confidence intervals.
- SciPy fitters use positive amplitudes, analytic Jacobians, two deterministic
  starts, tight tolerances, and at most five components. TRF uses bounded
  log rates/amplitudes; LM uses their unconstrained logarithmic parameterization.
  Known-K variants receive only the true count, never true rates/amplitudes.
- BIC variants evaluate K=1..5, stopping early only when reaching the common
  `1e-28*signal_energy` floor makes larger models incapable of improving BIC.
  The original residual-stop fitter is retained for comparability, not treated
  as the only SciPy alternative. Known-K oracle outcomes are reported below.
- Inverse references are analytic formulas evaluated at 60 digits then rounded
  to float64. Zero displayed error means agreement after rounding. mpmath's
  internal precision can exceed float64, so speed is not equal-precision proof.
- Forward callable tolerances match at 1e-10 absolute and 1e-8 relative.
  LaplaceX truncates at 1e6; omitted tails are negligible here. SciPy also gets
  an explicitly labeled squared-time endpoint variant for diffusion.
- These benchmarks use decaying real exponential mixtures, not arbitrary
  oscillatory signals or relaxation distributions. Decomposition no longer
  fabricates oscillation frequencies from Mexican-Hat side lobes. Manual
  oscillatory component reconstruction remains available.
- White homoscedastic noise is assumed. Colored noise, drift, unobserved
  time scales, and nearly coincident rates can invalidate count/rate inference.
  Fit reconstruction accuracy alone is not proof of uniquely identified modes.
- Configuration compatibility: `model_selection='bic'` is the new default;
  choose `'aicc'` for AICc. The old `use_aic` flag is retained as the enable
  switch for information-criterion pruning. `disable_als` now skips joint LM
  refinement despite its historical name. Defaults for scales/grid also changed.

## All forward measurements

| Case | Method | Median ms [min, max] | Max abs error | Target passed |
|---|---|---:|---:|---|
| exp_four | LaplaceX | 0.557 [0.528, 0.579] | 7.08e-16 | Yes |
| exp_four | SciPy quad | 0.946 [0.844, 1.163] | 1.39e-16 | Yes |
| exp_four | SciPy quad_vec | 4.692 [4.188, 4.794] | 5.55e-17 | Yes |
| exp_batch | LaplaceX | 8.871 [7.811, 9.744] | 2.05e-15 | Yes |
| exp_batch | SciPy quad | 5.740 [4.688, 8.686] | 2.65e-13 | Yes |
| exp_batch | SciPy quad_vec | 5.654 [4.767, 8.298] | 1.11e-16 | Yes |
| diffusion_four | LaplaceX | 0.962 [0.847, 1.021] | 6.90e-15 | Yes |
| diffusion_four | SciPy quad | 2.101 [2.067, 3.741] | 1.72e-12 | Yes |
| diffusion_four | SciPy quad_vec | Failed: ZeroDivisionError: division by zero | — | No |
| diffusion_four | SciPy quad_vec t=u^2 | 7.090 [6.883, 8.088] | 1.00e-15 | Yes |
| sampled_exp_2000 | LaplaceX spline | 7.138 [6.872, 7.674] | 2.41e-08 | No |
| sampled_exp_2000 | SciPy simpson | 1.321 [1.279, 1.589] | 6.64e-07 | No |
| sampled_exp_2000 | NumPy trapezoid | 0.775 [0.715, 0.975] | 1.90e-04 | No |
| sampled_exp_10000 | LaplaceX spline | 38.138 [36.752, 39.491] | 1.92e-10 | Yes |
| sampled_exp_10000 | SciPy simpson | 4.776 [4.531, 5.932] | 1.06e-09 | No |
| sampled_exp_10000 | NumPy trapezoid | 3.877 [3.442, 4.197] | 7.61e-06 | No |

## All inverse measurements

| Case | Method | Median ms [min, max] | Max abs error | Target passed |
|---|---|---:|---:|---|
| exp | LaplaceX N64 | 1.541 [1.328, 1.596] | 1.61e-14 | Yes |
| exp | LaplaceX N256 | 5.225 [5.014, 5.737] | 1.17e-14 | Yes |
| exp | mpmath Cohen degree32 | 80.790 [78.111, 82.336] | 2.74e-10 | No |
| exp | mpmath Cohen degree64 | 184.883 [177.809, 185.878] | 1.73e-18 | Yes |
| exp | mpmath Talbot degree32 | 279.604 [269.670, 284.865] | 8.67e-19 | Yes |
| exp | mpmath deHoog degree32 | 3416.472 [3278.216, 3679.011] | 0.00e+00 | Yes |
| exp | mpmath Stehfest degree32 | 317.723 [311.848, 329.291] | 6.04e-10 | No |
| exp | mpmath Cohen dps30 auto | 159.585 [158.388, 179.734] | 0.00e+00 | Yes |
| exp | SciPy impulse (rational coefficients) | 0.520 [0.493, 0.784] | 2.22e-16 | Yes |
| exp | LaplaceX N512 | 15.952 [10.762, 16.846] | 4.96e-14 | Yes |
| repeated_pole | LaplaceX N64 | 1.710 [1.407, 2.203] | 3.16e-15 | Yes |
| repeated_pole | LaplaceX N256 | 4.645 [4.379, 4.816] | 8.10e-15 | Yes |
| repeated_pole | mpmath Cohen degree32 | 112.054 [94.695, 151.140] | 3.76e-10 | Yes |
| repeated_pole | mpmath Cohen degree64 | 209.532 [201.895, 214.932] | 0.00e+00 | Yes |
| repeated_pole | mpmath Talbot degree32 | 282.473 [278.349, 330.538] | 0.00e+00 | Yes |
| repeated_pole | mpmath deHoog degree32 | 3346.282 [3264.145, 3553.061] | 0.00e+00 | Yes |
| repeated_pole | mpmath Stehfest degree32 | 332.719 [318.009, 420.661] | 3.78e-10 | Yes |
| repeated_pole | mpmath Cohen dps30 auto | 187.699 [185.381, 191.571] | 0.00e+00 | Yes |
| repeated_pole | SciPy impulse (rational coefficients) | 0.573 [0.528, 1.053] | 1.11e-16 | Yes |
| repeated_pole | LaplaceX N512 | 8.346 [6.584, 9.374] | 2.13e-14 | Yes |
| damped_cos | LaplaceX N64 | 4.777 [4.337, 6.450] | 3.32e-12 | Yes |
| damped_cos | LaplaceX N256 | 12.833 [11.630, 14.198] | 3.73e-14 | Yes |
| damped_cos | mpmath Cohen degree32 | 136.319 [126.560, 174.440] | 1.36e-06 | No |
| damped_cos | mpmath Cohen degree64 | 284.396 [275.574, 314.236] | 4.39e-15 | Yes |
| damped_cos | mpmath Talbot degree32 | 321.732 [313.221, 347.953] | 3.90e-09 | No |
| damped_cos | mpmath deHoog degree32 | 3729.587 [3332.968, 4612.636] | 0.00e+00 | Yes |
| damped_cos | mpmath Stehfest degree32 | 363.662 [350.082, 470.101] | 2.52e-02 | No |
| damped_cos | mpmath Cohen dps30 auto | 270.644 [245.783, 330.601] | 0.00e+00 | Yes |
| damped_cos | SciPy impulse (rational coefficients) | 0.544 [0.510, 0.752] | 4.47e-16 | Yes |
| damped_cos | LaplaceX N512 | 32.615 [26.803, 42.380] | 3.44e-14 | Yes |
| sin10 | LaplaceX N64 | 3.364 [2.907, 4.029] | 6.16e-02 | No |
| sin10 | LaplaceX N256 | 24.238 [23.687, 29.470] | 3.54e-09 | No |
| sin10 | mpmath Cohen degree32 | 97.777 [93.066, 110.911] | 5.86e+01 | No |
| sin10 | mpmath Cohen degree64 | 217.852 [216.172, 395.114] | 9.78e-06 | No |
| sin10 | mpmath Talbot degree32 | 296.873 [291.506, 355.178] | 1.51e+00 | No |
| sin10 | mpmath deHoog degree32 | 3419.225 [3268.961, 3435.064] | 1.45e-02 | No |
| sin10 | mpmath Stehfest degree32 | 324.447 [319.285, 425.793] | 1.03e+00 | No |
| sin10 | mpmath Cohen dps30 auto | 194.418 [183.933, 197.558] | 1.11e-16 | Yes |
| sin10 | SciPy impulse (rational coefficients) | 0.590 [0.579, 1.200] | 1.09e-14 | Yes |
| sin10 | LaplaceX N512 | 91.557 [80.770, 93.821] | 1.50e-14 | Yes |
| diffusion | LaplaceX N64 | 3.918 [2.776, 5.205] | 6.41e-01 | No |
| diffusion | LaplaceX N256 | 23.771 [22.737, 26.278] | 3.61e-01 | No |
| diffusion | mpmath Cohen degree32 | 134.371 [130.578, 152.530] | 1.29e-10 | Yes |
| diffusion | mpmath Cohen degree64 | 325.745 [294.188, 339.255] | 0.00e+00 | Yes |
| diffusion | mpmath Talbot degree32 | 360.084 [320.848, 421.573] | 0.00e+00 | Yes |
| diffusion | mpmath deHoog degree32 | 3886.281 [3719.212, 4418.348] | 0.00e+00 | Yes |
| diffusion | mpmath Stehfest degree32 | 352.137 [325.157, 432.129] | 6.66e-16 | Yes |
| diffusion | mpmath Cohen dps30 auto | 277.175 [268.602, 358.064] | 0.00e+00 | Yes |
| diffusion | LaplaceX N512 | 84.036 [83.111, 92.503] | 2.70e-01 | No |
| delayed_step | LaplaceX N64 | 4.260 [3.608, 4.978] | 1.35e-01 | No |
| delayed_step | LaplaceX N256 | 27.974 [26.202, 28.739] | 7.64e-02 | No |
| delayed_step | mpmath Cohen degree32 | 147.785 [142.615, 164.214] | 6.59e+02 | No |
| delayed_step | mpmath Cohen degree64 | 402.238 [386.401, 462.935] | 2.05e+07 | No |
| delayed_step | mpmath Talbot degree32 | Failed: ValueError: Output shape mismatch or nonfinite output | — | No |
| delayed_step | mpmath deHoog degree32 | 3742.342 [3509.587, 4373.065] | 6.89e-08 | No |
| delayed_step | mpmath Stehfest degree32 | 373.048 [350.305, 395.762] | 3.24e-01 | No |
| delayed_step | mpmath Cohen dps30 auto | 307.656 [302.546, 376.685] | 3.85e+00 | No |
| delayed_step | LaplaceX N512 | 92.712 [87.855, 106.501] | 3.07e-02 | No |
| delayed_step | mpmath deHoog degree64 | 15644.638 [15086.297, 18211.294] | 1.07e-13 | Yes |

## All component measurements

RMSE is against the clean generating signal. Rate error is the largest
relative error after sorting rates; it is omitted when counts disagree.

| Case | Method | Median ms [min, max] | Max abs error | Clean RMSE | Modes | Rate error |
|---|---|---:|---:|---:|---:|---:|
| single_uniform | LaplaceX default | 0.494 [0.355, 0.986] | 3.33e-16 | 3.23e-17 | 1 | 0.00e+00 |
| single_uniform | LaplaceX fine wavelets | 0.494 [0.342, 0.862] | 3.33e-16 | 3.23e-17 | 1 | 0.00e+00 |
| single_uniform | LaplaceX fine wavelets only | 0.208 [0.182, 0.375] | 3.33e-16 | 3.23e-17 | 1 | 0.00e+00 |
| single_uniform | LaplaceX no wavelets | 0.284 [0.212, 0.319] | 3.33e-16 | 3.23e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy least_squares auto K | 5.401 [5.115, 8.366] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy least_squares known K | 5.547 [5.012, 8.561] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy TRF auto K BIC | 5.701 [4.783, 7.206] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy LM auto K BIC | 3.034 [2.088, 4.358] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy LM known K | 2.425 [2.271, 3.611] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| two_uniform | LaplaceX default | 4.026 [3.912, 5.683] | 4.44e-16 | 9.98e-17 | 2 | 6.66e-16 |
| two_uniform | LaplaceX fine wavelets | 7.848 [7.260, 8.299] | 4.44e-16 | 9.98e-17 | 2 | 6.66e-16 |
| two_uniform | LaplaceX fine wavelets only | 7.334 [6.996, 8.296] | 4.44e-16 | 9.98e-17 | 2 | 6.66e-16 |
| two_uniform | LaplaceX no wavelets | 3.632 [3.556, 4.456] | 4.44e-16 | 4.53e-17 | 2 | 0.00e+00 |
| two_uniform | SciPy least_squares auto K | 12.622 [12.095, 15.713] | 4.44e-16 | 6.47e-17 | 2 | 1.11e-16 |
| two_uniform | SciPy least_squares known K | 6.744 [6.268, 8.899] | 4.44e-16 | 6.47e-17 | 2 | 1.11e-16 |
| two_uniform | SciPy TRF auto K BIC | 12.452 [11.497, 14.237] | 4.44e-16 | 6.47e-17 | 2 | 1.11e-16 |
| two_uniform | SciPy LM auto K BIC | 12.726 [11.389, 13.592] | 4.44e-16 | 6.72e-17 | 2 | 5.55e-16 |
| two_uniform | SciPy LM known K | 9.767 [9.184, 16.326] | 4.44e-16 | 6.72e-17 | 2 | 5.55e-16 |
| three_log_pilot | LaplaceX default | 2.873 [2.601, 4.068] | 9.99e-16 | 4.41e-16 | 3 | 2.66e-15 |
| three_log_pilot | LaplaceX fine wavelets | 4.384 [4.225, 5.476] | 9.99e-16 | 4.41e-16 | 3 | 2.66e-15 |
| three_log_pilot | LaplaceX fine wavelets only | 4.780 [3.783, 5.163] | 9.99e-16 | 4.41e-16 | 3 | 2.66e-15 |
| three_log_pilot | LaplaceX no wavelets | 25.600 [24.111, 28.240] | 4.44e-16 | 7.92e-17 | 3 | 1.11e-16 |
| three_log_pilot | SciPy least_squares auto K | 103.259 [87.870, 105.838] | 1.33e-15 | 3.08e-16 | 3 | 1.33e-15 |
| three_log_pilot | SciPy least_squares known K | 6.844 [6.380, 14.879] | 1.33e-15 | 3.08e-16 | 3 | 1.33e-15 |
| three_log_pilot | SciPy TRF auto K BIC | 98.538 [83.600, 107.759] | 1.33e-15 | 3.08e-16 | 3 | 1.33e-15 |
| three_log_pilot | SciPy LM auto K BIC | 37.752 [35.681, 40.273] | 6.66e-16 | 2.18e-16 | 3 | 4.44e-16 |
| three_log_pilot | SciPy LM known K | 4.748 [4.190, 6.730] | 6.66e-16 | 2.18e-16 | 3 | 4.44e-16 |
| two_log_pilot | LaplaceX default | 2.412 [2.292, 2.996] | 2.22e-16 | 2.77e-17 | 2 | 0.00e+00 |
| two_log_pilot | LaplaceX fine wavelets | 3.969 [3.769, 4.560] | 2.22e-16 | 2.77e-17 | 2 | 0.00e+00 |
| two_log_pilot | LaplaceX fine wavelets only | 3.670 [3.583, 4.550] | 2.22e-16 | 2.77e-17 | 2 | 0.00e+00 |
| two_log_pilot | LaplaceX no wavelets | 3.834 [3.492, 4.918] | 4.44e-16 | 8.47e-17 | 2 | 2.22e-16 |
| two_log_pilot | SciPy least_squares auto K | 12.483 [11.840, 12.960] | 4.44e-16 | 1.22e-16 | 2 | 5.55e-16 |
| two_log_pilot | SciPy least_squares known K | 5.575 [4.736, 5.621] | 4.44e-16 | 1.22e-16 | 2 | 5.55e-16 |
| two_log_pilot | SciPy TRF auto K BIC | 12.547 [11.465, 13.406] | 4.44e-16 | 1.22e-16 | 2 | 5.55e-16 |
| two_log_pilot | SciPy LM auto K BIC | 7.311 [6.540, 7.719] | 4.44e-16 | 1.24e-16 | 2 | 4.44e-16 |
| two_log_pilot | SciPy LM known K | 4.440 [4.228, 5.939] | 4.44e-16 | 1.24e-16 | 2 | 4.44e-16 |
| three_log_holdout | LaplaceX default | 2.863 [2.816, 3.742] | 4.44e-16 | 1.30e-16 | 3 | 6.66e-16 |
| three_log_holdout | LaplaceX fine wavelets | 4.906 [4.529, 5.723] | 4.44e-16 | 1.30e-16 | 3 | 6.66e-16 |
| three_log_holdout | LaplaceX fine wavelets only | 4.446 [4.346, 4.660] | 4.44e-16 | 1.30e-16 | 3 | 6.66e-16 |
| three_log_holdout | LaplaceX no wavelets | 10.372 [9.870, 11.255] | 1.11e-15 | 3.15e-16 | 3 | 5.55e-15 |
| three_log_holdout | SciPy least_squares auto K | 36.955 [35.959, 44.714] | 8.88e-16 | 1.72e-16 | 3 | 1.22e-15 |
| three_log_holdout | SciPy least_squares known K | 7.309 [7.069, 9.015] | 8.88e-16 | 1.72e-16 | 3 | 1.22e-15 |
| three_log_holdout | SciPy TRF auto K BIC | 38.353 [36.698, 45.609] | 8.88e-16 | 1.72e-16 | 3 | 1.22e-15 |
| three_log_holdout | SciPy LM auto K BIC | 17.693 [16.349, 18.268] | 4.44e-16 | 1.64e-16 | 3 | 2.22e-16 |
| three_log_holdout | SciPy LM known K | 5.007 [4.716, 5.313] | 4.44e-16 | 1.64e-16 | 3 | 2.22e-16 |
| close_rates_holdout | LaplaceX default | 12.872 [12.592, 13.099] | 4.44e-16 | 9.48e-17 | 2 | 5.87e-14 |
| close_rates_holdout | LaplaceX fine wavelets | 16.378 [15.629, 16.444] | 4.44e-16 | 9.48e-17 | 2 | 5.87e-14 |
| close_rates_holdout | LaplaceX fine wavelets only | 15.870 [15.163, 16.262] | 4.44e-16 | 9.48e-17 | 2 | 5.87e-14 |
| close_rates_holdout | LaplaceX no wavelets | 8.683 [8.560, 9.375] | 4.44e-16 | 7.12e-17 | 2 | 8.88e-16 |
| close_rates_holdout | SciPy least_squares auto K | 58.238 [57.045, 60.871] | 4.44e-16 | 6.45e-17 | 2 | 5.77e-15 |
| close_rates_holdout | SciPy least_squares known K | 52.939 [50.647, 58.461] | 4.44e-16 | 6.45e-17 | 2 | 5.77e-15 |
| close_rates_holdout | SciPy TRF auto K BIC | 57.513 [56.879, 65.300] | 4.44e-16 | 6.45e-17 | 2 | 5.77e-15 |
| close_rates_holdout | SciPy LM auto K BIC | 54.216 [50.675, 54.705] | 4.44e-16 | 5.81e-17 | 2 | 3.77e-15 |
| close_rates_holdout | SciPy LM known K | 49.251 [48.749, 50.991] | 4.44e-16 | 5.81e-17 | 2 | 3.77e-15 |
| weak_fast | LaplaceX default | 4.247 [4.042, 5.180] | 4.44e-16 | 1.51e-16 | 3 | 8.39e-13 |
| weak_fast | LaplaceX fine wavelets | 6.627 [6.001, 7.541] | 4.44e-16 | 1.51e-16 | 3 | 8.39e-13 |
| weak_fast | LaplaceX fine wavelets only | 5.662 [5.399, 6.012] | 4.44e-16 | 1.51e-16 | 3 | 8.39e-13 |
| weak_fast | LaplaceX no wavelets | 12.777 [11.207, 13.626] | 4.44e-16 | 1.82e-16 | 3 | 1.11e-14 |
| weak_fast | SciPy least_squares auto K | 27.518 [26.218, 33.899] | 4.25e-04 | 2.00e-04 | 2 | Count mismatch |
| weak_fast | SciPy least_squares known K | 16.469 [13.185, 29.832] | 4.44e-16 | 1.42e-16 | 3 | 6.62e-14 |
| weak_fast | SciPy TRF auto K BIC | 45.339 [39.792, 47.352] | 4.44e-16 | 1.42e-16 | 3 | 6.62e-14 |
| weak_fast | SciPy LM auto K BIC | 18.882 [16.370, 29.667] | 6.66e-16 | 1.42e-16 | 3 | 9.81e-14 |
| weak_fast | SciPy LM known K | 6.264 [5.736, 6.620] | 6.66e-16 | 1.42e-16 | 3 | 9.81e-14 |
| weak_middle | LaplaceX default | 3.759 [3.484, 4.428] | 4.44e-16 | 8.55e-17 | 3 | 2.69e-14 |
| weak_middle | LaplaceX fine wavelets | 5.393 [5.189, 5.822] | 4.44e-16 | 8.55e-17 | 3 | 2.69e-14 |
| weak_middle | LaplaceX fine wavelets only | 5.423 [4.753, 6.223] | 4.44e-16 | 8.55e-17 | 3 | 2.69e-14 |
| weak_middle | LaplaceX no wavelets | 10.082 [9.739, 11.038] | 4.44e-16 | 9.18e-17 | 3 | 6.99e-14 |
| weak_middle | SciPy least_squares auto K | 15.133 [14.701, 15.919] | 3.79e-04 | 1.93e-04 | 2 | Count mismatch |
| weak_middle | SciPy least_squares known K | 14.646 [13.918, 15.798] | 4.44e-16 | 1.20e-16 | 3 | 3.05e-14 |
| weak_middle | SciPy TRF auto K BIC | 29.707 [29.444, 31.807] | 4.44e-16 | 1.20e-16 | 3 | 3.05e-14 |
| weak_middle | SciPy LM auto K BIC | 13.573 [13.100, 14.060] | 4.44e-16 | 1.27e-16 | 3 | 2.05e-14 |
| weak_middle | SciPy LM known K | 6.499 [5.583, 6.819] | 4.44e-16 | 1.27e-16 | 3 | 2.05e-14 |
| weak_slow | LaplaceX default | 4.320 [3.429, 4.688] | 5.20e-16 | 2.24e-16 | 3 | 1.46e-12 |
| weak_slow | LaplaceX fine wavelets | 5.762 [5.125, 6.370] | 5.20e-16 | 2.24e-16 | 3 | 1.46e-12 |
| weak_slow | LaplaceX fine wavelets only | 5.575 [5.350, 6.472] | 5.20e-16 | 2.24e-16 | 3 | 1.46e-12 |
| weak_slow | LaplaceX no wavelets | 11.377 [10.249, 12.136] | 4.44e-16 | 1.35e-16 | 3 | 3.55e-15 |
| weak_slow | SciPy least_squares auto K | 27.511 [26.334, 28.804] | 8.53e-04 | 3.37e-04 | 2 | Count mismatch |
| weak_slow | SciPy least_squares known K | 12.982 [11.993, 13.417] | 4.44e-16 | 1.03e-16 | 3 | 1.54e-14 |
| weak_slow | SciPy TRF auto K BIC | 42.018 [39.990, 49.222] | 4.44e-16 | 1.03e-16 | 3 | 1.54e-14 |
| weak_slow | SciPy LM auto K BIC | 86.749 [82.942, 97.726] | 4.44e-16 | 9.89e-17 | 3 | 1.11e-15 |
| weak_slow | SciPy LM known K | 65.105 [63.738, 78.310] | 4.44e-16 | 9.89e-17 | 3 | 1.11e-15 |
| noisy_log_seed7 | LaplaceX default | 4.762 [4.366, 4.864] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | LaplaceX fine wavelets | 7.168 [6.136, 8.580] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | LaplaceX fine wavelets only | 6.964 [6.436, 7.918] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | LaplaceX no wavelets | 12.224 [12.143, 13.079] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy least_squares auto K | 39.316 [36.954, 47.103] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy least_squares known K | 9.519 [8.550, 10.510] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy TRF auto K BIC | 550.669 [529.327, 762.510] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy LM auto K BIC | 459.299 [448.855, 564.905] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy LM known K | 5.274 [4.446, 5.707] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed19 | LaplaceX default | 4.382 [4.094, 4.845] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | LaplaceX fine wavelets | 6.640 [6.205, 7.090] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | LaplaceX fine wavelets only | 6.150 [5.653, 7.454] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | LaplaceX no wavelets | 13.951 [12.916, 14.538] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy least_squares auto K | 37.599 [36.177, 42.427] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy least_squares known K | 6.988 [6.738, 7.777] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy TRF auto K BIC | 831.837 [795.172, 856.437] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy LM auto K BIC | 553.416 [544.316, 599.431] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy LM known K | 5.110 [4.786, 5.413] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed41 | LaplaceX default | 4.846 [4.692, 5.468] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | LaplaceX fine wavelets | 7.005 [6.630, 7.561] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | LaplaceX fine wavelets only | 6.563 [6.371, 7.778] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | LaplaceX no wavelets | 13.345 [12.412, 13.849] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy least_squares auto K | 36.930 [35.755, 50.283] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy least_squares known K | 7.346 [6.853, 12.010] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy TRF auto K BIC | 548.284 [543.757, 589.059] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy LM auto K BIC | 442.327 [440.263, 458.834] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy LM known K | 4.917 [4.884, 5.747] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |

## Reproduce

```bash
python -m pip install ".[benchmark]"
python benchmarks/validate_wavelets.py
python benchmarks/compare.py --section all --repeats 5 --points 64
python benchmarks/compare_extended.py --repeats 5 --points 64
python benchmarks/render_comparison.py --update-readme
```

Data: [forward](comparison/forward.json), [inverse](comparison/inverse.json),
[inverse high-order follow-up](comparison/inverse_extended.json),
[components](comparison/decomposition.json), [validation](comparison/validation.json).

## Primary references

- [MINPACK QR/Levenberg–Marquardt algorithms](https://www.netlib.org/minpack/)
- [SciPy least_squares: TRF and LM](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)
- [mpmath inversion algorithms and precision](https://mpmath.org/doc/current/calculus/inverselaplace.html)
- [SciPy quad](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html)
- [SciPy quad_vec](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad_vec.html)
- [SciPy impulse](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.impulse.html)
