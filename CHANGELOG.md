# Changelog

All notable changes to LaplaceX are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

---

## [1.2.0] — 2026-10-21

- Refine rates by separable least squares with QR-fitted amplitudes and an exact
  projected Jacobian. Retain the original joint optimizer as a retry when a
  projected result violates the existing coalescence/amplitude limits.
- Search unresolved pairs even without an above-noise residual ridge. Use a
  permissive residual/noise screening bound before the existing BIC/AICc test.
- Add detection regression fixtures, an independently seeded noise/separation
  sweep, signed/null checks, and matched before/after timing controls.
- Update `benchmarks/COMPARISON.md` with accuracy, speed, and remaining failures.

## [1.1.0] — 2026-10-21

- Add `ExponentialComponents.weights` for normalized absolute amplitudes;
  retain `significances` as a compatibility alias. These are not statistical
  significance measures. The numerical core is unchanged.
- Add external compiled log/linear-coordinate controls, variable-projection
  and matrix-pencil comparisons, a frozen noise/separation sweep with local
  Gaussian Cramer--Rao diagnostics, and attributed experimental fluorescence data.
- Add independent checks for the evaluation Jacobian, Fisher information,
  specialist baselines, and amplitude-weight compatibility.

## [1.0.0] — 2026-10-21

### Noisy component detection
- Enforce conditional amplitude identifiability for noisy fast modes whose
  characteristic time precedes the observed window; refit the retained modes.
- Recalibrate automatic proposal noise with residual MAD once the fit reaches
  the noise level, recovering weak modes hidden by the initial estimate.
- Add exact-observation regressions for both failures and a real fast-mode
  preservation test; validate 480 mixtures across three fixed seeds.

### Pair-case performance follow-up
- Cache propagated wavelet noise variance only at candidate crossings and pad
  convolution boundaries once per residual search.
- Refine unresolved-pair proposals on both sides of their residual ridge using
  conditional Jacobian projections, retaining the original proposals.
- Reuse a Householder reduction across LM damping trials; both wavelet and
  rate-grid fits share this solver optimization without relaxed tolerances.
- Add eight clean/noisy signed-pair regressions, rerun both validation seeds,
  and preserve the previous measurements for the targeted before/after comparison.

### Log-plane wavelet rewrite
- Track persistent Mexican-Hat zero crossings from fine to coarse scales on
  `f(exp(x))`, with exact propagation of interpolation/filter noise covariance.
- Initialize persistent modes jointly and replace coordinate rate searches with
  an analytic-Jacobian QR/Levenberg–Marquardt solver on the original samples.
- Add BIC selection, conditional amplitude uncertainty and refitted deletions;
  expose `model_selection="bic"` or `"aicc"`. `use_aic` remains the pruning enable switch.
- Change default initialization to 16 scales over [0.05, 0.5] and 256 log points.
- Enforce mode budgets in both stages, handle constants and noise-only signals,
  and avoid redundant noise-driven fallback fits and residual-energy overflow.
- Restrict noise-covariance accumulation to the filter's sample support and share
  a numerically verified single-exponential shortcut between both fitting routes.
- Add weak-mode, scaling, noise and mode-count regressions, a separate 160-case
  validation, stronger SciPy LM/BIC comparisons, and preserved historical results.

### Fixed
- Corrected the embedded Gauss weights, restored integration at zero, enforced
  `t_max`, and made quadrature fail explicitly when it cannot converge.
- Replaced the mislabeled Stehfest inverse with actual Weeks–Laguerre inversion;
  `sigma`, `b`, and `N` now control the stated method.
- Corrected conjugate poles and phase in oscillatory component transforms.
- Fixed NumPy slice handling, preserved evaluation shapes, and validated sample
  ordering, lengths, finite values, and numerical configuration.
- Explicit inverse callables take precedence over cached decompositions.
- Fixed AIC scaling and exact-fit scores; component significances are normalized
  after pruning. Unsupported wavelets raise an error.
- Release wheel tests now run the suite and propagate failures.

### Improved
- Precomputed wavelet kernels; reused fitting basis columns; replaced normal
  equations with Householder QR and converged rate refinement.
- Sampled transforms integrate each spline segment, extrapolate the first
  segment to zero, and use zero continuation after the final sample.
- Added analytic accuracy regressions and reproducible before/after benchmarks.
- Portable builds no longer enable `-march=native` by default; opt in with
  `-DLAPLACEX_NATIVE=ON` for machine-specific builds.

## Initial development release — 2026-05-23

### Added
- **C++ core** (`laplacex::core`) — adaptive Gauss-Kronrod G7/K15 quadrature
  for the forward transform with exponential axis compression; parallelised
  batch evaluation via `std::async`.
- **Cubic-spline forward transform** from discrete samples.
- **Weeks–Laguerre inverse** with automatic parameter selection.
- **Critical-points module** (`laplacex::critical`) — continuous wavelet
  transform (Mexican Hat) in the logarithmic time plane; ridge-line extraction;
  AIC-based component pruning.
- **Reconstruction module** — time-domain and Laplace-domain reconstruction
  from `ExponentialComponent` lists.
- **Utilities** — test-signal generators, noise estimation, log-plane
  resampling, AIC/BIC criteria.
- **pybind11 bindings** — full C++ surface exposed to Python with NumPy
  array support.
- **Python API** — `laplace()`, `inverse()`, `decompose()`, `reconstruct()`,
  `LaplaceTransform`, `ExponentialComponents`.
- **CI/CD** — GitHub Actions matrix (Linux, macOS Intel, macOS ARM, Windows),
  cibuildwheel for pre-built wheels, OIDC PyPI publishing.
- **Full pytest suite** — 40+ tests covering correctness, edge cases, and
  round-trip behaviour.
