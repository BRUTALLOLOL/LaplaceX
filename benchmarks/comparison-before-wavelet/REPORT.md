# Comparison with Python alternatives

Measured locally on 20 October 2026. Every published time is the median of
5 measured runs after one warmup; raw run times, failures,
warnings, input definitions, versions, and hashes are preserved in JSON.

![Comparison and wavelet ablation](comparison/comparison.png)

## Findings

LaplaceX has a workload-specific speed advantage over mpmath for batches
of smooth real inverse-transform values. The advantage comes from a compiled
Weeks expansion whose coefficients are reused across times. That inverse
path does not perform wavelet analysis. This is not evidence that wavelets
are responsible for the inverse-transform speed advantage.

| Inverse workload (64 times) | LaplaceX, ms / max error | Fastest passing mpmath, ms / max error | Speed ratio |
|---|---:|---:|---:|
| exp | 1.29 / 1.61e-14 (LaplaceX N64) | 157.12 / 0.00e+00 (mpmath Cohen dps30 auto) | 121.4× |
| repeated_pole | 1.46 / 3.16e-15 (LaplaceX N64) | 93.56 / 3.76e-10 (mpmath Cohen degree32) | 64.2× |
| damped_cos | 4.41 / 3.32e-12 (LaplaceX N64) | 258.78 / 0.00e+00 (mpmath Cohen dps30 auto) | 58.7× |
| sin10 | 87.15 / 1.50e-14 (LaplaceX N512) | 192.56 / 1.11e-16 (mpmath Cohen dps30 auto) | 2.2× |
| diffusion | No tested setting passed | 133.03 / 1.29e-10 (mpmath Cohen degree32) | — |
| delayed_step | No tested setting passed | 14830.83 / 1.07e-13 (mpmath deHoog degree64) | — |

The speed ratio compares the fastest **passing** tested configuration in
each package, selected separately for each case. Passing means every point
satisfies `abs(error) <= 1e-10 + 1e-8*abs(reference)`. No speed ratio is
claimed if either package has no passing configuration. Full results below
include all tested configurations and failures, including difficult cases.

SciPy is competitive or faster on several forward and fitting workloads.
A general claim of speed or accuracy supremacy is not supported.

## What the wavelet ablation establishes

The ablation bypasses CWT by calling the existing native `decompose_deflate`
entry point directly. It keeps the same compiled QR/ALS fitting code, sample
arrays, residual stopping tolerance, and AIC setting. Fallback fitting is
limited to five components; CWT determines its own initial mode count.
The no-wavelet route starts its own rate grid search; it receives no true
rates or amplitudes. This tests the value of wavelet initialization relative
to the existing fallback, not an unrelated implementation language. It removes
the log-CWT stage as a whole; it does not separately isolate logarithmic
coordinates from the wavelet operation.

Default CWT often finds only one component on mixtures; the public API
then repeats the no-wavelet fallback. Identical results with extra CWT
cost are evidence against attributing those fits to wavelet analysis.

The finer preset (`scale_min=0.05`, `scale_max=0.5`, `n_scales=32`,
`min_ridge_len=3`) was explored on the `two_uniform`, `three_log_pilot`,
and `two_log_pilot` signals before the final run. It is an exploratory
configuration, not a universally tuned default. Additional rates, a close
pair, and three noise seeds test its transfer beyond those pilots.

On some three-component log-grid cases the fine-wavelet route avoids the
independent rate-search fallback and modestly reduces fitting time. The
counterexamples matter: single- and two-component cases can be slower or
over-segmented, and noisy-fit errors can be effectively identical with or
without wavelets. SciPy fitting remains faster on the tested mixtures.

Logarithmic wavelets are therefore a configurable initialization mechanism
with local measured gains, not a demonstrated general accuracy advantage.

## Fairness and precision

- Python: `3.14.7 (tags/v3.14.7:823f032, Aug  5 2026, 10:51:32) [MSC v.1944 64 bit (AMD64)]`; platform: `Windows-10-10.0.19045-SP0`.
- NumPy `2.5.3`, SciPy `1.18.1`, mpmath `1.3.0`.
- Benchmark sections run serially; each uses one BLAS thread. Method order
  is shuffled deterministically on each repetition. Setup and array
  conversion inside each method are timed; imports and truth generation are not.
- LaplaceX uses the compiled Release wheel from the preceding fixes:
  Zig 0.16.0 / Clang 21.1.0, C++17, `-O3`, Windows GNU x64 target.
  The extension SHA-256 and benchmark-script hash are recorded in each JSON.
- Inverse references use analytic formulas evaluated at 60 decimal digits,
  then rounded to float64. All outputs are compared as float64. This
  measures float64 workloads; mpmath can deliver much higher precision.
  A displayed zero error means agreement after float64 rounding, not an
  exact numerical inversion.
- mpmath includes Cohen at degrees 32 and 64, Talbot/de Hoog/Stehfest at
  degree 32, and automatic Cohen at requested 30 dps. Fixed-degree calls
  begin at 15 dps but mpmath raises its working precision internally.
  Equal term counts do not imply equal work or equal arithmetic precision.
  The local mpmath backend is `python`; optional gmpy2 acceleration
  was not installed. Results do not compare every possible package build.
- A separate follow-up adds LaplaceX N512 on every inverse case and de Hoog
  degree64 on the delayed step after the initial settings exposed failures.
  Original runs are retained; follow-up settings are labeled and saved separately.
- Both LaplaceX and mpmath inverse calls receive only a callable image.
  `scipy.signal.impulse` receives numerator/denominator coefficients and
  is explicitly labeled as a richer-input rational-system comparison.
- Forward callable comparisons use shared absolute/relative tolerances
  `1e-10`/`1e-8`. SciPy integrates to infinity; LaplaceX truncates at `1e6`.
  Omitted tails are negligible for these specific frequencies and signals.
  The diffusion test also gives SciPy the same squared-time endpoint
  regularization used by LaplaceX, as an explicitly labeled variant.
- Sampled-data methods receive the same samples and finite-window truth.
  They are evaluated for speed and approximation error, not treated as
  interchangeable algorithms with identical accuracy guarantees.
- Decomposition judges predictions against the clean signal even for
  noisy inputs; it also reports true/recovered counts and rate errors.
  SciPy fitting uses positive amplitudes, an analytic Jacobian, bounded
  log rates, two deterministic starts, and tight optimizer tolerances.
  Automatic K uses the same residual-energy stopping rule and five-mode
  budget. Known-K receives the true count only, not rates or amplitudes.
- These synthetic cases and three noise seeds are a bounded local study,
  not confidence intervals for all inputs, CPUs, or real experimental data.

## All forward measurements

| Case | Method | Median ms [min, max] | Max abs error | Target passed |
|---|---|---:|---:|---|
| exp_four | LaplaceX | 0.567 [0.476, 0.866] | 7.08e-16 | Yes |
| exp_four | SciPy quad | 0.802 [0.749, 1.089] | 1.39e-16 | Yes |
| exp_four | SciPy quad_vec | 4.280 [4.007, 5.012] | 5.55e-17 | Yes |
| exp_batch | LaplaceX | 7.436 [7.147, 8.317] | 2.05e-15 | Yes |
| exp_batch | SciPy quad | 4.786 [4.583, 5.354] | 2.65e-13 | Yes |
| exp_batch | SciPy quad_vec | 4.954 [4.767, 5.666] | 1.11e-16 | Yes |
| diffusion_four | LaplaceX | 0.877 [0.814, 0.911] | 6.90e-15 | Yes |
| diffusion_four | SciPy quad | 2.280 [2.141, 2.789] | 1.72e-12 | Yes |
| diffusion_four | SciPy quad_vec | Failed: ZeroDivisionError: division by zero | — | No |
| diffusion_four | SciPy quad_vec t=u^2 | 6.714 [6.588, 9.061] | 1.00e-15 | Yes |
| sampled_exp_2000 | LaplaceX spline | 7.126 [6.669, 7.491] | 2.41e-08 | No |
| sampled_exp_2000 | SciPy simpson | 1.209 [1.121, 1.986] | 6.64e-07 | No |
| sampled_exp_2000 | NumPy trapezoid | 0.655 [0.652, 1.006] | 1.90e-04 | No |
| sampled_exp_10000 | LaplaceX spline | 37.529 [36.906, 39.193] | 1.92e-10 | Yes |
| sampled_exp_10000 | SciPy simpson | 4.628 [4.442, 5.474] | 1.06e-09 | No |
| sampled_exp_10000 | NumPy trapezoid | 3.604 [3.402, 3.974] | 7.61e-06 | No |

## All inverse measurements

| Case | Method | Median ms [min, max] | Max abs error | Target passed |
|---|---|---:|---:|---|
| exp | LaplaceX N64 | 1.294 [1.261, 1.429] | 1.61e-14 | Yes |
| exp | LaplaceX N256 | 4.989 [4.954, 5.653] | 1.17e-14 | Yes |
| exp | mpmath Cohen degree32 | 82.265 [77.314, 89.532] | 2.74e-10 | No |
| exp | mpmath Cohen degree64 | 179.677 [173.755, 184.701] | 1.73e-18 | Yes |
| exp | mpmath Talbot degree32 | 264.225 [262.803, 267.657] | 8.67e-19 | Yes |
| exp | mpmath deHoog degree32 | 3230.272 [3176.909, 3340.386] | 0.00e+00 | Yes |
| exp | mpmath Stehfest degree32 | 312.262 [303.219, 328.090] | 6.04e-10 | No |
| exp | mpmath Cohen dps30 auto | 157.122 [149.106, 159.791] | 0.00e+00 | Yes |
| exp | SciPy impulse (rational coefficients) | 0.494 [0.481, 0.711] | 2.22e-16 | Yes |
| exp | LaplaceX N512 | 10.474 [9.678, 10.798] | 4.96e-14 | Yes |
| repeated_pole | LaplaceX N64 | 1.458 [1.411, 1.843] | 3.16e-15 | Yes |
| repeated_pole | LaplaceX N256 | 4.916 [4.581, 5.535] | 8.10e-15 | Yes |
| repeated_pole | mpmath Cohen degree32 | 93.561 [92.872, 102.301] | 3.76e-10 | Yes |
| repeated_pole | mpmath Cohen degree64 | 212.035 [204.794, 237.243] | 0.00e+00 | Yes |
| repeated_pole | mpmath Talbot degree32 | 280.743 [277.090, 301.916] | 0.00e+00 | Yes |
| repeated_pole | mpmath deHoog degree32 | 3311.875 [3224.265, 3646.628] | 0.00e+00 | Yes |
| repeated_pole | mpmath Stehfest degree32 | 323.469 [300.930, 335.245] | 3.78e-10 | Yes |
| repeated_pole | mpmath Cohen dps30 auto | 187.423 [179.669, 193.284] | 0.00e+00 | Yes |
| repeated_pole | SciPy impulse (rational coefficients) | 0.558 [0.532, 0.576] | 1.11e-16 | Yes |
| repeated_pole | LaplaceX N512 | 5.563 [5.261, 5.759] | 2.13e-14 | Yes |
| damped_cos | LaplaceX N64 | 4.405 [4.062, 4.431] | 3.32e-12 | Yes |
| damped_cos | LaplaceX N256 | 12.336 [11.238, 13.726] | 3.73e-14 | Yes |
| damped_cos | mpmath Cohen degree32 | 126.220 [124.216, 152.783] | 1.36e-06 | No |
| damped_cos | mpmath Cohen degree64 | 280.137 [277.640, 297.085] | 4.39e-15 | Yes |
| damped_cos | mpmath Talbot degree32 | 333.675 [306.177, 358.787] | 3.90e-09 | No |
| damped_cos | mpmath deHoog degree32 | 3365.373 [3273.600, 3746.114] | 0.00e+00 | Yes |
| damped_cos | mpmath Stehfest degree32 | 347.775 [342.163, 426.253] | 2.52e-02 | No |
| damped_cos | mpmath Cohen dps30 auto | 258.776 [238.286, 290.055] | 0.00e+00 | Yes |
| damped_cos | SciPy impulse (rational coefficients) | 0.514 [0.509, 0.598] | 4.47e-16 | Yes |
| damped_cos | LaplaceX N512 | 28.463 [27.248, 29.285] | 3.44e-14 | Yes |
| sin10 | LaplaceX N64 | 3.120 [3.034, 3.881] | 6.16e-02 | No |
| sin10 | LaplaceX N256 | 24.548 [23.418, 25.574] | 3.54e-09 | No |
| sin10 | mpmath Cohen degree32 | 94.818 [92.333, 131.050] | 5.86e+01 | No |
| sin10 | mpmath Cohen degree64 | 215.651 [212.006, 231.874] | 9.78e-06 | No |
| sin10 | mpmath Talbot degree32 | 292.493 [287.187, 314.497] | 1.51e+00 | No |
| sin10 | mpmath deHoog degree32 | 3326.158 [3291.440, 3370.172] | 1.45e-02 | No |
| sin10 | mpmath Stehfest degree32 | 318.393 [313.535, 394.019] | 1.03e+00 | No |
| sin10 | mpmath Cohen dps30 auto | 192.563 [183.761, 202.419] | 1.11e-16 | Yes |
| sin10 | SciPy impulse (rational coefficients) | 0.676 [0.579, 0.998] | 1.09e-14 | Yes |
| sin10 | LaplaceX N512 | 87.151 [85.923, 91.908] | 1.50e-14 | Yes |
| diffusion | LaplaceX N64 | 2.900 [2.695, 7.138] | 6.41e-01 | No |
| diffusion | LaplaceX N256 | 22.885 [22.477, 24.213] | 3.61e-01 | No |
| diffusion | mpmath Cohen degree32 | 133.033 [119.565, 146.179] | 1.29e-10 | Yes |
| diffusion | mpmath Cohen degree64 | 341.664 [291.991, 380.267] | 0.00e+00 | Yes |
| diffusion | mpmath Talbot degree32 | 317.230 [306.417, 410.113] | 0.00e+00 | Yes |
| diffusion | mpmath deHoog degree32 | 3349.940 [3281.628, 3453.340] | 0.00e+00 | Yes |
| diffusion | mpmath Stehfest degree32 | 332.191 [317.139, 392.962] | 6.66e-16 | Yes |
| diffusion | mpmath Cohen dps30 auto | 261.933 [249.020, 300.186] | 0.00e+00 | Yes |
| diffusion | LaplaceX N512 | 81.678 [73.080, 89.318] | 2.70e-01 | No |
| delayed_step | LaplaceX N64 | 3.754 [3.615, 5.262] | 1.35e-01 | No |
| delayed_step | LaplaceX N256 | 27.340 [25.884, 29.020] | 7.64e-02 | No |
| delayed_step | mpmath Cohen degree32 | 151.490 [138.633, 166.826] | 6.59e+02 | No |
| delayed_step | mpmath Cohen degree64 | 386.453 [372.420, 411.476] | 2.05e+07 | No |
| delayed_step | mpmath Talbot degree32 | Failed: ValueError: Output shape mismatch or nonfinite output | — | No |
| delayed_step | mpmath deHoog degree32 | 3392.378 [3331.347, 3471.289] | 6.89e-08 | No |
| delayed_step | mpmath Stehfest degree32 | 335.156 [324.558, 356.804] | 3.24e-01 | No |
| delayed_step | mpmath Cohen dps30 auto | 313.042 [303.366, 322.997] | 3.85e+00 | No |
| delayed_step | LaplaceX N512 | 89.436 [83.462, 93.267] | 3.07e-02 | No |
| delayed_step | mpmath deHoog degree64 | 14830.833 [14747.934, 15589.490] | 1.07e-13 | Yes |

## All decomposition measurements

Rate error is the largest relative difference after sorting true and fitted
rates; it is omitted when component counts disagree. Clean RMSE is computed
against the generating signal. `pilot` labels identify explored cases.

| Case | Method | Median ms [min, max] | Max abs error | Clean RMSE | Modes | Rate error |
|---|---|---:|---:|---:|---:|---:|
| single_uniform | LaplaceX default | 25.315 [24.747, 25.672] | 7.26e-14 | 1.07e-14 | 1 | 1.10e-13 |
| single_uniform | LaplaceX fine wavelets | 38.187 [36.923, 39.444] | 4.88e-15 | 3.62e-16 | 3 | Count mismatch |
| single_uniform | LaplaceX fine wavelets only | 37.973 [37.047, 39.657] | 4.88e-15 | 3.62e-16 | 3 | Count mismatch |
| single_uniform | LaplaceX no wavelets | 7.098 [6.311, 7.524] | 2.00e-14 | 2.23e-15 | 1 | 2.22e-14 |
| single_uniform | SciPy least_squares auto K | 5.005 [4.620, 5.427] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| single_uniform | SciPy least_squares known K | 4.945 [4.870, 6.117] | 2.22e-16 | 2.51e-17 | 1 | 0.00e+00 |
| two_uniform | LaplaceX default | 294.209 [289.811, 308.126] | 3.36e-12 | 1.31e-12 | 2 | 4.28e-11 |
| two_uniform | LaplaceX fine wavelets | 260.544 [260.031, 269.743] | 2.29e-12 | 8.81e-13 | 2 | 2.89e-11 |
| two_uniform | LaplaceX fine wavelets only | 261.819 [254.932, 266.017] | 2.29e-12 | 8.81e-13 | 2 | 2.89e-11 |
| two_uniform | LaplaceX no wavelets | 269.217 [263.530, 284.335] | 3.36e-12 | 1.31e-12 | 2 | 4.28e-11 |
| two_uniform | SciPy least_squares auto K | 12.261 [11.456, 13.563] | 4.44e-16 | 6.47e-17 | 2 | 1.11e-16 |
| two_uniform | SciPy least_squares known K | 6.493 [6.356, 8.758] | 4.44e-16 | 6.47e-17 | 2 | 1.11e-16 |
| three_log_pilot | LaplaceX default | 628.574 [619.482, 652.286] | 2.54e-12 | 9.28e-13 | 3 | 8.73e-12 |
| three_log_pilot | LaplaceX fine wavelets | 501.239 [487.089, 505.294] | 1.27e-12 | 4.90e-13 | 3 | 5.83e-12 |
| three_log_pilot | LaplaceX fine wavelets only | 497.097 [482.588, 505.895] | 1.27e-12 | 4.90e-13 | 3 | 5.83e-12 |
| three_log_pilot | LaplaceX no wavelets | 590.083 [580.543, 615.869] | 2.54e-12 | 9.28e-13 | 3 | 8.73e-12 |
| three_log_pilot | SciPy least_squares auto K | 90.811 [82.383, 95.010] | 1.33e-15 | 3.08e-16 | 3 | 1.33e-15 |
| three_log_pilot | SciPy least_squares known K | 6.410 [5.710, 10.715] | 1.33e-15 | 3.08e-16 | 3 | 1.33e-15 |
| two_log_pilot | LaplaceX default | 113.836 [110.370, 114.838] | 1.19e-13 | 4.37e-14 | 2 | 3.50e-13 |
| two_log_pilot | LaplaceX fine wavelets | 271.924 [265.505, 272.816] | 7.46e-13 | 2.34e-13 | 3 | Count mismatch |
| two_log_pilot | LaplaceX fine wavelets only | 270.671 [269.532, 276.143] | 7.46e-13 | 2.34e-13 | 3 | Count mismatch |
| two_log_pilot | LaplaceX no wavelets | 85.977 [84.564, 87.863] | 1.19e-13 | 4.37e-14 | 2 | 3.50e-13 |
| two_log_pilot | SciPy least_squares auto K | 12.387 [11.979, 15.757] | 4.44e-16 | 1.22e-16 | 2 | 5.55e-16 |
| two_log_pilot | SciPy least_squares known K | 5.142 [4.844, 5.955] | 4.44e-16 | 1.22e-16 | 2 | 5.55e-16 |
| three_log_holdout | LaplaceX default | 646.351 [638.643, 684.834] | 5.26e-12 | 2.07e-12 | 3 | 3.28e-11 |
| three_log_holdout | LaplaceX fine wavelets | 580.924 [579.972, 608.350] | 2.44e-12 | 9.59e-13 | 3 | 1.51e-11 |
| three_log_holdout | LaplaceX fine wavelets only | 579.559 [570.124, 657.821] | 2.44e-12 | 9.59e-13 | 3 | 1.51e-11 |
| three_log_holdout | LaplaceX no wavelets | 616.739 [613.259, 622.994] | 5.26e-12 | 2.07e-12 | 3 | 3.28e-11 |
| three_log_holdout | SciPy least_squares auto K | 37.247 [35.953, 37.876] | 8.88e-16 | 1.72e-16 | 3 | 1.22e-15 |
| three_log_holdout | SciPy least_squares known K | 7.068 [7.040, 8.213] | 8.88e-16 | 1.72e-16 | 3 | 1.22e-15 |
| close_rates_holdout | LaplaceX default | 1004.653 [997.678, 1053.697] | 8.62e-05 | 2.09e-05 | 2 | 1.59e-02 |
| close_rates_holdout | LaplaceX fine wavelets | 2001.167 [1951.522, 2065.344] | 6.04e-06 | 1.47e-06 | 3 | Count mismatch |
| close_rates_holdout | LaplaceX fine wavelets only | 2002.058 [1975.048, 2043.447] | 6.04e-06 | 1.47e-06 | 3 | Count mismatch |
| close_rates_holdout | LaplaceX no wavelets | 969.164 [964.351, 987.870] | 8.62e-05 | 2.09e-05 | 2 | 1.59e-02 |
| close_rates_holdout | SciPy least_squares auto K | 58.710 [56.505, 59.910] | 4.44e-16 | 6.45e-17 | 2 | 5.77e-15 |
| close_rates_holdout | SciPy least_squares known K | 52.804 [50.553, 57.370] | 4.44e-16 | 6.45e-17 | 2 | 5.77e-15 |
| noisy_log_seed7 | LaplaceX default | 574.989 [570.565, 587.462] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | LaplaceX fine wavelets | 513.181 [512.238, 523.251] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | LaplaceX fine wavelets only | 513.239 [501.198, 534.551] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | LaplaceX no wavelets | 544.020 [532.870, 552.308] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy least_squares auto K | 37.646 [35.283, 38.612] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed7 | SciPy least_squares known K | 8.201 [7.885, 12.556] | 3.74e-04 | 1.69e-04 | 3 | 1.65e-03 |
| noisy_log_seed19 | LaplaceX default | 581.832 [571.105, 599.963] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | LaplaceX fine wavelets | 498.673 [481.217, 522.208] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | LaplaceX fine wavelets only | 486.581 [486.309, 503.447] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | LaplaceX no wavelets | 544.718 [536.454, 565.015] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy least_squares auto K | 36.770 [34.845, 42.695] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed19 | SciPy least_squares known K | 7.444 [6.956, 8.470] | 1.58e-04 | 6.64e-05 | 3 | 8.41e-04 |
| noisy_log_seed41 | LaplaceX default | 582.785 [576.485, 621.155] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | LaplaceX fine wavelets | 513.125 [506.664, 517.943] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | LaplaceX fine wavelets only | 514.426 [507.148, 526.728] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | LaplaceX no wavelets | 539.280 [533.735, 555.106] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy least_squares auto K | 38.598 [36.533, 40.864] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |
| noisy_log_seed41 | SciPy least_squares known K | 7.217 [6.933, 12.885] | 1.77e-04 | 6.78e-05 | 3 | 1.63e-03 |

## Reproduce

```bash
python -m pip install ".[benchmark]"
python benchmarks/compare.py --section all --repeats 5 --points 64
python benchmarks/compare_extended.py --repeats 5 --points 64
python benchmarks/render_comparison.py
```

Raw data: [forward](comparison/forward.json), [inverse](comparison/inverse.json),
[inverse follow-up](comparison/inverse_extended.json),
[decomposition](comparison/decomposition.json). The script does not change
library code or automatically select favorable cases. Existing output files
are replaced when the corresponding section is rerun.

## Primary documentation

- [mpmath inverse algorithms and internal precision](https://mpmath.org/doc/current/calculus/inverselaplace.html)
- [SciPy quad tolerances and complex integrands](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html)
- [SciPy quad_vec](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad_vec.html)
- [SciPy least_squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)
- [SciPy rational-system impulse response](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.impulse.html)
