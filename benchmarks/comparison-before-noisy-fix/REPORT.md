# Log-plane wavelet recovery and Python comparisons

Measured on 21 October 2026 after the pair-case optimization. Times are medians
of five runs after one warmup, with shuffled method order and one BLAS thread.

![Wavelet recovery comparison](comparison/wavelet_comparison.png)

## Component recovery

| Signal | LaplaceX default, ms | Fastest passing SciPy auto-K, ms | Speed ratio | Modes: LaplaceX / truth |
|---|---:|---:|---:|---:|
| single_uniform | 0.41 | 2.20 | 5.37× | 1 / 1 |
| two_uniform | 2.49 | 11.00 | 4.41× | 2 / 2 |
| three_log_pilot | 1.18 | 35.16 | 29.90× | 3 / 3 |
| two_log_pilot | 0.98 | 6.98 | 7.15× | 2 / 2 |
| three_log_holdout | 1.29 | 17.67 | 13.67× | 3 / 3 |
| close_rates_holdout | 7.06 | 55.82 | 7.90× | 2 / 2 |
| weak_fast | 2.32 | 17.24 | 7.44× | 3 / 3 |
| weak_middle | 2.14 | 15.34 | 7.17× | 3 / 3 |
| weak_slow | 2.50 | 42.72 | 17.06× | 3 / 3 |
| noisy_log_seed7 | 2.63 | 35.53 | 13.52× | 3 / 3 |
| noisy_log_seed19 | 2.47 | 33.64 | 13.61× | 3 / 3 |
| noisy_log_seed41 | 2.50 | 37.12 | 14.84× | 3 / 3 |

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

## Optimizing the previously slower pair cases

- Evaluate and cache exact propagated noise variance only at candidate crossings.
- Pad convolution boundaries once per residual rather than clamp each filter tap.
- For unresolved pairs, search both sides of the strongest residual ridge using
  conditional Jacobian projections. Keep the original proposals and existing selection rules.
- Reduce the Jacobian with Householder QR once per nonlinear iteration, then
  try damping changes on the small triangular problem. No normal equations or
  relaxed tolerances are introduced. The QR improvement is shared by both routes.

| Previously slower case | Before: wavelet / grid ms | Now: wavelet / grid ms | Wavelet speedup vs previous | Now: grid / wavelet |
|---|---:|---:|---:|---:|
| two_uniform | 4.03 / 3.63 | 2.49 / 3.89 | 1.61× | 1.56× |
| close_rates_holdout | 12.87 / 8.68 | 7.06 / 8.63 | 1.82× | 1.22× |

A separate **21-run** randomized check of the final binary confirms the local pair timings.
Cells show median [25th, 75th percentile] milliseconds; every fit must pass the component-count and accuracy checks.

| Case | Default wavelets, ms | Rate grid, ms | Grid / wavelet |
|---|---:|---:|---:|
| two_uniform | 2.50 [2.45, 2.63] | 3.75 [3.60, 4.10] | 1.50× |
| close_rates_holdout | 6.81 [6.37, 7.25] | 7.73 [7.54, 8.70] | 1.14× |

Before/after timings come from separate local runs; same-run grid/wavelet ratios
are the cleaner ablation. All models and final fits still use the original observations.

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
Maximum clean rate error: `1.828e-12` relative.
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

A further seed (20261107) recovered **158/160** counts;
the pre-optimization algorithm recovered **158/160** on those same inputs.
There were **0 changed count decisions** on that fresh set. The two noisy failures
are preserved in the raw results; this optimization does not solve all count-identifiability limits.
The fresh set's noisy rate-error 95th percentile is **4.92% before / 4.92% after**.

The fresh-seed baseline was rebuilt from the archived
[previous critical.cpp](comparison-before-pair-optimization/critical.cpp) with the
same Release toolchain and Python wrapper; it was not inferred from the new fits.

## Wavelet ablation

`LaplaceX no wavelets` bypasses CWT and searches a 64-point log-rate grid.
It shares the new QR/LM solver, noise estimation, conditional candidate test,
BIC and deletion refits. This isolates initialization and candidate generation
from the optimizer rewrite. It removes the entire log-CWT stage; it does not
separately compare linear-time and log-time wavelets.

The default public route also includes Python fallback/noise checks. The
native no-wavelet entry omits that routing overhead; the explicitly
wavelet-only rows omit fallback checks and provide an additional comparison.

With the same optimizer, default wavelet initialization is **22.07×** faster
on the widely separated three-mode case and **7.65×** faster on its holdout.
The two-uniform and close-rate ratios are **1.56×** and
**1.22×** respectively; above 1 means wavelets are faster.
These ratios compare the native rate-grid route
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
| single_uniform | 25.32 / 1 | 0.41 / 1 | 7.26e-14 → 3.33e-16 |
| two_uniform | 294.21 / 2 | 2.49 / 2 | 3.36e-12 → 6.66e-16 |
| three_log_pilot | 628.57 / 3 | 1.18 / 3 | 2.54e-12 → 9.99e-16 |
| two_log_pilot | 113.84 / 2 | 0.98 / 2 | 1.19e-13 → 2.22e-16 |
| three_log_holdout | 646.35 / 3 | 1.29 / 3 | 5.26e-12 → 4.44e-16 |
| close_rates_holdout | 1004.65 / 2 | 7.06 / 2 | 8.62e-05 → 4.44e-16 |
| noisy_log_seed7 | 574.99 / 3 | 2.63 / 3 | 3.74e-04 → 3.74e-04 |
| noisy_log_seed19 | 581.83 / 3 | 2.47 / 3 | 1.58e-04 → 1.58e-04 |
| noisy_log_seed41 | 582.79 / 3 | 2.50 / 3 | 1.77e-04 → 1.77e-04 |

Historical raw files and benchmark scripts are preserved in
[comparison-before-wavelet](comparison-before-wavelet/decomposition.json).

## Forward and inverse comparisons

| Inverse workload (64 times) | LaplaceX, ms / max error | Fastest passing mpmath, ms / max error | Speed ratio |
|---|---:|---:|---:|
| exp | 1.39 / 1.61e-14 (LaplaceX N64) | 172.43 / 1.73e-18 (mpmath Cohen degree64) | 124.3× |
| repeated_pole | 1.49 / 3.16e-15 (LaplaceX N64) | 94.07 / 3.76e-10 (mpmath Cohen degree32) | 63.2× |
| damped_cos | 3.99 / 3.32e-12 (LaplaceX N64) | 259.48 / 0.00e+00 (mpmath Cohen dps30 auto) | 65.0× |
| sin10 | 98.07 / 1.50e-14 (LaplaceX N512) | 177.40 / 1.11e-16 (mpmath Cohen dps30 auto) | 1.8× |
| diffusion | No tested setting passed | 130.83 / 1.29e-10 (mpmath Cohen degree32) | — |
| delayed_step | No tested setting passed | 14371.24 / 1.07e-13 (mpmath deHoog degree64) | — |

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
| exp_four | LaplaceX | 0.603 [0.566, 1.744] | 7.08e-16 | Yes |
| exp_four | SciPy quad | 0.908 [0.896, 1.955] | 1.39e-16 | Yes |
| exp_four | SciPy quad_vec | 4.945 [4.319, 7.232] | 5.55e-17 | Yes |
| exp_batch | LaplaceX | 10.437 [9.567, 11.779] | 2.05e-15 | Yes |
| exp_batch | SciPy quad | 6.923 [4.552, 9.626] | 2.65e-13 | Yes |
| exp_batch | SciPy quad_vec | 7.472 [5.654, 9.746] | 1.11e-16 | Yes |
| diffusion_four | LaplaceX | 1.006 [0.918, 1.143] | 6.90e-15 | Yes |
| diffusion_four | SciPy quad | 3.415 [2.151, 5.882] | 1.72e-12 | Yes |
| diffusion_four | SciPy quad_vec | Failed: ZeroDivisionError: division by zero | — | No |
| diffusion_four | SciPy quad_vec t=u^2 | 8.861 [6.878, 10.428] | 1.00e-15 | Yes |
| sampled_exp_2000 | LaplaceX spline | 7.906 [7.121, 8.612] | 2.41e-08 | No |
| sampled_exp_2000 | SciPy simpson | 1.377 [1.269, 1.926] | 6.64e-07 | No |
| sampled_exp_2000 | NumPy trapezoid | 0.882 [0.689, 1.076] | 1.90e-04 | No |
| sampled_exp_10000 | LaplaceX spline | 42.314 [39.183, 44.327] | 1.92e-10 | Yes |
| sampled_exp_10000 | SciPy simpson | 5.982 [4.989, 7.232] | 1.06e-09 | No |
| sampled_exp_10000 | NumPy trapezoid | 4.310 [3.965, 5.112] | 7.61e-06 | No |

## All inverse measurements

| Case | Method | Median ms [min, max] | Max abs error | Target passed |
|---|---|---:|---:|---|
| exp | LaplaceX N64 | 1.387 [1.258, 1.691] | 1.61e-14 | Yes |
| exp | LaplaceX N256 | 5.464 [5.154, 7.241] | 1.17e-14 | Yes |
| exp | mpmath Cohen degree32 | 89.560 [75.115, 96.704] | 2.74e-10 | No |
| exp | mpmath Cohen degree64 | 172.428 [165.880, 201.398] | 1.73e-18 | Yes |
| exp | mpmath Talbot degree32 | 286.522 [246.980, 343.451] | 8.67e-19 | Yes |
| exp | mpmath deHoog degree32 | 3501.959 [3290.300, 3965.062] | 0.00e+00 | Yes |
| exp | mpmath Stehfest degree32 | 328.304 [284.540, 416.082] | 6.04e-10 | No |
| exp | mpmath Cohen dps30 auto | 175.188 [142.697, 233.309] | 0.00e+00 | Yes |
| exp | SciPy impulse (rational coefficients) | 0.512 [0.479, 0.770] | 2.22e-16 | Yes |
| exp | LaplaceX N512 | 9.531 [8.872, 12.793] | 4.96e-14 | Yes |
| repeated_pole | LaplaceX N64 | 1.488 [1.382, 2.470] | 3.16e-15 | Yes |
| repeated_pole | LaplaceX N256 | 4.591 [4.367, 9.120] | 8.10e-15 | Yes |
| repeated_pole | mpmath Cohen degree32 | 94.071 [83.192, 206.652] | 3.76e-10 | Yes |
| repeated_pole | mpmath Cohen degree64 | 232.546 [192.495, 313.411] | 0.00e+00 | Yes |
| repeated_pole | mpmath Talbot degree32 | 320.373 [272.620, 488.465] | 0.00e+00 | Yes |
| repeated_pole | mpmath deHoog degree32 | 3416.960 [2991.468, 4275.197] | 0.00e+00 | Yes |
| repeated_pole | mpmath Stehfest degree32 | 363.304 [288.363, 436.075] | 3.78e-10 | Yes |
| repeated_pole | mpmath Cohen dps30 auto | 208.200 [172.287, 296.157] | 0.00e+00 | Yes |
| repeated_pole | SciPy impulse (rational coefficients) | 0.526 [0.521, 0.558] | 1.11e-16 | Yes |
| repeated_pole | LaplaceX N512 | 5.887 [5.460, 6.200] | 2.13e-14 | Yes |
| damped_cos | LaplaceX N64 | 3.990 [3.924, 5.780] | 3.32e-12 | Yes |
| damped_cos | LaplaceX N256 | 10.869 [10.693, 14.491] | 3.73e-14 | Yes |
| damped_cos | mpmath Cohen degree32 | 116.693 [110.541, 122.316] | 1.36e-06 | No |
| damped_cos | mpmath Cohen degree64 | 262.203 [248.931, 375.030] | 4.39e-15 | Yes |
| damped_cos | mpmath Talbot degree32 | 303.612 [287.606, 334.654] | 3.90e-09 | No |
| damped_cos | mpmath deHoog degree32 | 3132.265 [3081.853, 3691.117] | 0.00e+00 | Yes |
| damped_cos | mpmath Stehfest degree32 | 318.593 [305.053, 329.819] | 2.52e-02 | No |
| damped_cos | mpmath Cohen dps30 auto | 259.484 [231.602, 277.480] | 0.00e+00 | Yes |
| damped_cos | SciPy impulse (rational coefficients) | 0.519 [0.483, 0.945] | 4.47e-16 | Yes |
| damped_cos | LaplaceX N512 | 28.005 [25.935, 37.810] | 3.44e-14 | Yes |
| sin10 | LaplaceX N64 | 3.817 [2.979, 4.563] | 6.16e-02 | No |
| sin10 | LaplaceX N256 | 24.271 [23.117, 39.545] | 3.54e-09 | No |
| sin10 | mpmath Cohen degree32 | 98.170 [86.033, 140.465] | 5.86e+01 | No |
| sin10 | mpmath Cohen degree64 | 213.199 [191.607, 309.281] | 9.78e-06 | No |
| sin10 | mpmath Talbot degree32 | 285.592 [270.299, 350.320] | 1.51e+00 | No |
| sin10 | mpmath deHoog degree32 | 3367.519 [3234.094, 3435.273] | 1.45e-02 | No |
| sin10 | mpmath Stehfest degree32 | 367.147 [291.657, 464.914] | 1.03e+00 | No |
| sin10 | mpmath Cohen dps30 auto | 177.398 [172.378, 204.170] | 1.11e-16 | Yes |
| sin10 | SciPy impulse (rational coefficients) | 0.591 [0.528, 0.859] | 1.09e-14 | Yes |
| sin10 | LaplaceX N512 | 98.073 [79.014, 108.713] | 1.50e-14 | Yes |
| diffusion | LaplaceX N64 | 3.071 [2.702, 3.782] | 6.41e-01 | No |
| diffusion | LaplaceX N256 | 22.298 [21.888, 23.390] | 3.61e-01 | No |
| diffusion | mpmath Cohen degree32 | 130.829 [129.161, 133.427] | 1.29e-10 | Yes |
| diffusion | mpmath Cohen degree64 | 312.592 [285.052, 404.287] | 0.00e+00 | Yes |
| diffusion | mpmath Talbot degree32 | 336.263 [301.968, 367.785] | 0.00e+00 | Yes |
| diffusion | mpmath deHoog degree32 | 3487.591 [3356.733, 3975.456] | 0.00e+00 | Yes |
| diffusion | mpmath Stehfest degree32 | 343.278 [311.212, 374.228] | 6.66e-16 | Yes |
| diffusion | mpmath Cohen dps30 auto | 268.035 [254.998, 289.965] | 0.00e+00 | Yes |
| diffusion | LaplaceX N512 | 76.742 [73.011, 78.005] | 2.70e-01 | No |
| delayed_step | LaplaceX N64 | 4.679 [3.443, 5.472] | 1.35e-01 | No |
| delayed_step | LaplaceX N256 | 28.319 [26.204, 29.420] | 7.64e-02 | No |
| delayed_step | mpmath Cohen degree32 | 180.151 [163.831, 181.303] | 6.59e+02 | No |
| delayed_step | mpmath Cohen degree64 | 409.345 [383.545, 563.347] | 2.05e+07 | No |
| delayed_step | mpmath Talbot degree32 | Failed: ValueError: Output shape mismatch or nonfinite output | — | No |
| delayed_step | mpmath deHoog degree32 | 3540.834 [3274.323, 4081.812] | 6.89e-08 | No |
| delayed_step | mpmath Stehfest degree32 | 399.694 [327.285, 500.731] | 3.24e-01 | No |
| delayed_step | mpmath Cohen dps30 auto | 313.289 [287.623, 353.794] | 3.85e+00 | No |
| delayed_step | LaplaceX N512 | 90.805 [79.313, 96.302] | 3.07e-02 | No |
| delayed_step | mpmath deHoog degree64 | 14371.241 [12335.066, 17117.394] | 1.07e-13 | Yes |

## All component measurements

RMSE is against the clean generating signal. Rate error is the largest
relative error after sorting rates; it is omitted when counts disagree.

| Case | Method | Median ms [min, max] | Max abs error | Clean RMSE | Modes | Rate error |
|---|---|---:|---:|---:|---:|---:|
| single_uniform | LaplaceX default | 0.409 [0.398, 0.433] | 3.33e-16 | 3.23e-17 | 1 | 0.00e+00 |
| single_uniform | LaplaceX fine wavelets | 0.378 [0.349, 0.455] | 3.33e-16 | 3.23e-17 | 1 | 0.00e+00 |
| single_uniform | LaplaceX fine wavelets only | 0.204 [0.191, 0.414] | 3.33e-16 | 3.23e-17 | 1 | 0.00e+00 |
| single_uniform | LaplaceX no wavelets | 0.235 [0.181, 0.297] | 3.33e-16 | 3.23e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy least_squares auto K | 5.551 [5.043, 7.804] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy least_squares known K | 4.753 [4.644, 6.231] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy TRF auto K BIC | 5.722 [4.716, 7.543] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy LM auto K BIC | 2.197 [2.001, 4.432] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy LM known K | 2.093 [1.958, 2.252] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| two_uniform | LaplaceX default | 2.494 [2.345, 3.201] | 6.66e-16 | 1.37e-16 | 2 | 4.44e-16 |
| two_uniform | LaplaceX fine wavelets | 3.380 [3.022, 3.955] | 6.66e-16 | 1.37e-16 | 2 | 4.44e-16 |
| two_uniform | LaplaceX fine wavelets only | 3.236 [2.802, 3.645] | 6.66e-16 | 1.37e-16 | 2 | 4.44e-16 |
| two_uniform | LaplaceX no wavelets | 3.891 [3.649, 4.221] | 4.44e-16 | 4.05e-17 | 2 | 2.22e-16 |
| two_uniform | SciPy least_squares auto K | 11.078 [10.541, 13.622] | 4.44e-16 | 6.47e-17 | 2 | 1.11e-16 |
| two_uniform | SciPy least_squares known K | 7.505 [5.611, 8.846] | 4.44e-16 | 6.47e-17 | 2 | 1.11e-16 |
| two_uniform | SciPy TRF auto K BIC | 11.467 [10.610, 15.020] | 4.44e-16 | 6.47e-17 | 2 | 1.11e-16 |
| two_uniform | SciPy LM auto K BIC | 10.996 [10.159, 12.668] | 4.44e-16 | 6.72e-17 | 2 | 5.55e-16 |
| two_uniform | SciPy LM known K | 10.343 [9.180, 11.977] | 4.44e-16 | 6.72e-17 | 2 | 5.55e-16 |
| three_log_pilot | LaplaceX default | 1.176 [1.116, 1.600] | 9.99e-16 | 4.41e-16 | 3 | 2.66e-15 |
| three_log_pilot | LaplaceX fine wavelets | 1.608 [1.406, 1.890] | 9.99e-16 | 4.41e-16 | 3 | 2.66e-15 |
| three_log_pilot | LaplaceX fine wavelets only | 1.314 [1.121, 1.867] | 9.99e-16 | 4.41e-16 | 3 | 2.66e-15 |
| three_log_pilot | LaplaceX no wavelets | 25.958 [22.386, 27.176] | 4.44e-16 | 7.92e-17 | 3 | 1.11e-16 |
| three_log_pilot | SciPy least_squares auto K | 85.614 [78.906, 92.101] | 1.33e-15 | 3.08e-16 | 3 | 1.33e-15 |
| three_log_pilot | SciPy least_squares known K | 5.588 [5.347, 8.616] | 1.33e-15 | 3.08e-16 | 3 | 1.33e-15 |
| three_log_pilot | SciPy TRF auto K BIC | 89.530 [78.980, 98.756] | 1.33e-15 | 3.08e-16 | 3 | 1.33e-15 |
| three_log_pilot | SciPy LM auto K BIC | 35.163 [32.543, 51.238] | 6.66e-16 | 2.18e-16 | 3 | 4.44e-16 |
| three_log_pilot | SciPy LM known K | 4.767 [3.796, 5.315] | 6.66e-16 | 2.18e-16 | 3 | 4.44e-16 |
| two_log_pilot | LaplaceX default | 0.976 [0.896, 1.136] | 2.22e-16 | 2.77e-17 | 2 | 0.00e+00 |
| two_log_pilot | LaplaceX fine wavelets | 1.151 [1.133, 1.792] | 2.22e-16 | 2.77e-17 | 2 | 0.00e+00 |
| two_log_pilot | LaplaceX fine wavelets only | 0.940 [0.904, 1.014] | 2.22e-16 | 2.77e-17 | 2 | 0.00e+00 |
| two_log_pilot | LaplaceX no wavelets | 3.845 [3.796, 4.173] | 4.44e-16 | 8.47e-17 | 2 | 2.22e-16 |
| two_log_pilot | SciPy least_squares auto K | 12.138 [11.282, 12.934] | 4.44e-16 | 1.22e-16 | 2 | 5.55e-16 |
| two_log_pilot | SciPy least_squares known K | 5.546 [4.902, 6.102] | 4.44e-16 | 1.22e-16 | 2 | 5.55e-16 |
| two_log_pilot | SciPy TRF auto K BIC | 13.438 [12.290, 16.666] | 4.44e-16 | 1.22e-16 | 2 | 5.55e-16 |
| two_log_pilot | SciPy LM auto K BIC | 6.975 [6.412, 7.665] | 4.44e-16 | 1.24e-16 | 2 | 4.44e-16 |
| two_log_pilot | SciPy LM known K | 4.474 [4.395, 5.269] | 4.44e-16 | 1.24e-16 | 2 | 4.44e-16 |
| three_log_holdout | LaplaceX default | 1.293 [1.204, 1.522] | 4.44e-16 | 1.30e-16 | 3 | 6.66e-16 |
| three_log_holdout | LaplaceX fine wavelets | 1.465 [1.400, 1.678] | 4.44e-16 | 1.30e-16 | 3 | 6.66e-16 |
| three_log_holdout | LaplaceX fine wavelets only | 1.637 [1.244, 1.681] | 4.44e-16 | 1.30e-16 | 3 | 6.66e-16 |
| three_log_holdout | LaplaceX no wavelets | 9.886 [9.234, 10.786] | 1.11e-15 | 3.20e-16 | 3 | 5.55e-15 |
| three_log_holdout | SciPy least_squares auto K | 36.060 [33.589, 50.019] | 8.88e-16 | 1.72e-16 | 3 | 1.22e-15 |
| three_log_holdout | SciPy least_squares known K | 7.689 [6.403, 8.456] | 8.88e-16 | 1.72e-16 | 3 | 1.22e-15 |
| three_log_holdout | SciPy TRF auto K BIC | 36.936 [33.400, 37.720] | 8.88e-16 | 1.72e-16 | 3 | 1.22e-15 |
| three_log_holdout | SciPy LM auto K BIC | 17.674 [16.713, 19.030] | 4.44e-16 | 1.64e-16 | 3 | 2.22e-16 |
| three_log_holdout | SciPy LM known K | 4.808 [4.501, 5.589] | 4.44e-16 | 1.64e-16 | 3 | 2.22e-16 |
| close_rates_holdout | LaplaceX default | 7.062 [6.048, 8.171] | 4.44e-16 | 6.23e-17 | 2 | 1.11e-15 |
| close_rates_holdout | LaplaceX fine wavelets | 7.815 [7.113, 9.146] | 4.44e-16 | 6.23e-17 | 2 | 1.11e-15 |
| close_rates_holdout | LaplaceX fine wavelets only | 7.239 [6.958, 8.023] | 4.44e-16 | 6.23e-17 | 2 | 1.11e-15 |
| close_rates_holdout | LaplaceX no wavelets | 8.633 [7.956, 10.622] | 4.44e-16 | 6.99e-17 | 2 | 2.44e-15 |
| close_rates_holdout | SciPy least_squares auto K | 62.542 [55.724, 77.906] | 4.44e-16 | 6.45e-17 | 2 | 5.77e-15 |
| close_rates_holdout | SciPy least_squares known K | 56.802 [52.791, 63.278] | 4.44e-16 | 6.45e-17 | 2 | 5.77e-15 |
| close_rates_holdout | SciPy TRF auto K BIC | 60.661 [52.884, 72.198] | 4.44e-16 | 6.45e-17 | 2 | 5.77e-15 |
| close_rates_holdout | SciPy LM auto K BIC | 55.825 [50.738, 61.820] | 4.44e-16 | 5.81e-17 | 2 | 3.77e-15 |
| close_rates_holdout | SciPy LM known K | 52.249 [48.359, 57.354] | 4.44e-16 | 5.81e-17 | 2 | 3.77e-15 |
| weak_fast | LaplaceX default | 2.319 [2.252, 2.662] | 4.44e-16 | 1.52e-16 | 3 | 8.35e-13 |
| weak_fast | LaplaceX fine wavelets | 2.649 [2.576, 2.679] | 4.44e-16 | 1.52e-16 | 3 | 8.35e-13 |
| weak_fast | LaplaceX fine wavelets only | 2.485 [2.463, 2.556] | 4.44e-16 | 1.52e-16 | 3 | 8.35e-13 |
| weak_fast | LaplaceX no wavelets | 10.335 [9.945, 11.921] | 4.44e-16 | 1.82e-16 | 3 | 1.09e-14 |
| weak_fast | SciPy least_squares auto K | 26.040 [25.098, 29.268] | 4.25e-04 | 2.00e-04 | 2 | Count mismatch |
| weak_fast | SciPy least_squares known K | 13.644 [12.829, 14.199] | 4.44e-16 | 1.42e-16 | 3 | 6.62e-14 |
| weak_fast | SciPy TRF auto K BIC | 46.326 [38.477, 50.364] | 4.44e-16 | 1.42e-16 | 3 | 6.62e-14 |
| weak_fast | SciPy LM auto K BIC | 17.239 [16.584, 23.946] | 6.66e-16 | 1.42e-16 | 3 | 9.81e-14 |
| weak_fast | SciPy LM known K | 6.347 [6.010, 7.429] | 6.66e-16 | 1.42e-16 | 3 | 9.81e-14 |
| weak_middle | LaplaceX default | 2.140 [2.000, 2.527] | 4.44e-16 | 8.54e-17 | 3 | 3.46e-14 |
| weak_middle | LaplaceX fine wavelets | 2.713 [2.437, 3.604] | 4.44e-16 | 8.54e-17 | 3 | 3.46e-14 |
| weak_middle | LaplaceX fine wavelets only | 2.357 [2.232, 2.983] | 4.44e-16 | 8.54e-17 | 3 | 3.46e-14 |
| weak_middle | LaplaceX no wavelets | 8.932 [8.461, 9.033] | 4.44e-16 | 9.25e-17 | 3 | 8.79e-14 |
| weak_middle | SciPy least_squares auto K | 15.258 [14.762, 22.885] | 3.79e-04 | 1.93e-04 | 2 | Count mismatch |
| weak_middle | SciPy least_squares known K | 15.287 [13.891, 15.940] | 4.44e-16 | 1.20e-16 | 3 | 3.05e-14 |
| weak_middle | SciPy TRF auto K BIC | 30.690 [28.049, 35.132] | 4.44e-16 | 1.20e-16 | 3 | 3.05e-14 |
| weak_middle | SciPy LM auto K BIC | 15.343 [12.884, 20.386] | 4.44e-16 | 1.27e-16 | 3 | 2.05e-14 |
| weak_middle | SciPy LM known K | 6.070 [5.549, 7.462] | 4.44e-16 | 1.27e-16 | 3 | 2.05e-14 |
| weak_slow | LaplaceX default | 2.504 [2.236, 3.148] | 5.20e-16 | 2.24e-16 | 3 | 1.46e-12 |
| weak_slow | LaplaceX fine wavelets | 2.565 [2.517, 2.936] | 5.20e-16 | 2.24e-16 | 3 | 1.46e-12 |
| weak_slow | LaplaceX fine wavelets only | 2.591 [2.367, 3.094] | 5.20e-16 | 2.24e-16 | 3 | 1.46e-12 |
| weak_slow | LaplaceX no wavelets | 10.189 [9.866, 11.240] | 4.44e-16 | 1.35e-16 | 3 | 3.55e-15 |
| weak_slow | SciPy least_squares auto K | 28.467 [25.563, 42.252] | 8.53e-04 | 3.37e-04 | 2 | Count mismatch |
| weak_slow | SciPy least_squares known K | 12.652 [11.932, 13.239] | 4.44e-16 | 1.03e-16 | 3 | 1.54e-14 |
| weak_slow | SciPy TRF auto K BIC | 42.717 [38.168, 44.598] | 4.44e-16 | 1.03e-16 | 3 | 1.54e-14 |
| weak_slow | SciPy LM auto K BIC | 81.069 [79.295, 88.265] | 4.44e-16 | 9.89e-17 | 3 | 1.11e-15 |
| weak_slow | SciPy LM known K | 68.288 [62.833, 87.230] | 4.44e-16 | 9.89e-17 | 3 | 1.11e-15 |
| noisy_log_seed7 | LaplaceX default | 2.628 [2.526, 3.147] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | LaplaceX fine wavelets | 3.479 [3.406, 3.936] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | LaplaceX fine wavelets only | 2.999 [2.909, 3.155] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | LaplaceX no wavelets | 12.797 [11.763, 13.016] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy least_squares auto K | 35.526 [34.780, 39.141] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy least_squares known K | 7.696 [7.458, 8.761] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy TRF auto K BIC | 527.697 [512.348, 562.096] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy LM auto K BIC | 464.402 [420.928, 562.696] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy LM known K | 4.732 [4.628, 4.880] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed19 | LaplaceX default | 2.472 [2.367, 2.745] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | LaplaceX fine wavelets | 3.183 [2.969, 4.193] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | LaplaceX fine wavelets only | 2.825 [2.653, 3.234] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | LaplaceX no wavelets | 11.610 [11.457, 14.746] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy least_squares auto K | 33.637 [32.720, 35.689] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy least_squares known K | 6.785 [6.555, 7.122] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy TRF auto K BIC | 761.504 [758.265, 878.207] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy LM auto K BIC | 524.926 [502.302, 569.221] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy LM known K | 4.717 [4.391, 5.280] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed41 | LaplaceX default | 2.502 [2.392, 2.650] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | LaplaceX fine wavelets | 3.140 [2.868, 3.716] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | LaplaceX fine wavelets only | 2.841 [2.654, 3.366] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | LaplaceX no wavelets | 13.081 [11.667, 28.584] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy least_squares auto K | 37.124 [33.545, 39.867] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy least_squares known K | 7.265 [6.737, 8.298] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy TRF auto K BIC | 530.075 [502.014, 546.428] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy LM auto K BIC | 448.383 [410.200, 454.247] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy LM known K | 4.638 [4.313, 5.076] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |

## Reproduce

```bash
python -m pip install ".[benchmark]"
python benchmarks/validate_wavelets.py
python benchmarks/validate_wavelets.py --seed 20261107 --output benchmarks/comparison/validation_fresh.json
python benchmarks/compare.py --section all --repeats 5 --points 64
python benchmarks/compare_extended.py --repeats 5 --points 64
python benchmarks/compare_pairs.py --repeats 21
python benchmarks/render_comparison.py --update-readme
```

Data: [forward](comparison/forward.json), [inverse](comparison/inverse.json),
[inverse high-order follow-up](comparison/inverse_extended.json),
[components](comparison/decomposition.json), [validation](comparison/validation.json),
[fresh validation](comparison/validation_fresh.json),
[previous pair-case measurements](comparison-before-pair-optimization/decomposition.json),
[previous algorithm on the fresh seed](comparison-before-pair-optimization/validation_fresh.json).

[Longer pair timing check](comparison/pairs.json).

## Primary references

- [MINPACK QR/Levenberg–Marquardt algorithms](https://www.netlib.org/minpack/)
- [SciPy least_squares: TRF and LM](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)
- [mpmath inversion algorithms and precision](https://mpmath.org/doc/current/calculus/inverselaplace.html)
- [SciPy quad](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html)
- [SciPy quad_vec](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad_vec.html)
- [SciPy impulse](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.impulse.html)
