# LaplaceX 1.0.0

**Python library for numerical forward and inverse Laplace transforms.**

LaplaceX combines adaptive Gauss–Kronrod forward quadrature, Weeks–Laguerre
inversion, and exponential fitting seeded by wavelet analysis in logarithmic
time. Accuracy depends on the signal, sampling window, and configuration;
see [measured accuracy and limitations](benchmarks/ACCURACY.md).

| Operation | Implementation |
|---|---|
| Callable forward transform | Adaptive embedded G7–K15 quadrature |
| Sampled forward transform | Natural cubic spline, integrated at every knot |
| Inverse transform | Weeks–Laguerre expansion with automatic scale search |
| Exponential decomposition | Persistent log-wavelet ridges, joint QR/LM fitting, BIC |
| Python interface | NumPy arrays and a compiled C++ core via pybind11 |

---


**Requirements:** Python ≥ 3.10, NumPy ≥ 1.24.
Optional extras: `matplotlib` (plotting), `pandas` (DataFrame export).

---

<!-- comparison:start -->
## Log-wavelet recovery: measured results

Persistent log-plane wavelet ridges now initialize a joint rate/amplitude fit.
Noise propagation, BIC selection and refitted deletions control missed and
duplicate modes. The default route does not need the true component count.

The separately seeded validation recovered the correct count in **160/160** mixtures
with the rate-grid fallback disabled. It includes 80 clean and 80 noisy signals,
1–4 modes, linear and logarithmic sampling, and independently varied amplitudes/rates.
Maximum clean rate error: `1.828e-12` relative.
For noise σ=0.001, median / 95th-percentile / maximum rate errors were
**0.26% / 4.55% / 7.30%**.
Correct counts do not imply exact rates in noisy, poorly conditioned mixtures.
These are finite synthetic validation results, not a guarantee for all signals.

The two noisy failures in seed 20261107 are fixed: **158/160 → 160/160** correct counts.
There are **2 changed count decisions** on that regression set.
An additional seed (20261119) recovers **160/160** counts.
Across the three seeds: **480/480**, with fallback disabled.
The formerly failing set's noisy rate-error 95th percentile is **4.92% before / 4.85% after**.
The additional seed's noisy rate-error 95th percentile is **6.58%**.
Across all noisy mixtures, the largest rate error is **21.07%**;
this pre-existing limitation shows why correct counts do not guarantee precise individual rates.
Seed 20261107 is now regression data because its failures informed these fixes;
seed 20261119 was chosen before implementation and was not used to tune the rules.

The boundary uncertainty check removes an unconstrained extrapolated mode.
Residual MAD calibration reveals a previously suppressed weak component.
Explicit noise settings are respected, and identifiable fast modes remain supported.

| Previously failing case | Modes before → after / truth | Clean RMSE before → after | Largest rate error after |
|---|---:|---:|---:|
| log grid, case 41 | 3 → 2 / 2 | 1.78e-04 → 9.72e-05 | 0.17% |
| uniform grid, case 47 | 3 → 4 / 4 | 2.71e-04 → 9.28e-05 | 2.93% |

![Component speed and weak-mode accuracy](benchmarks/comparison/wavelet_comparison.png)

| Signal | LaplaceX default, ms | Fastest passing SciPy auto-K, ms | Speed ratio | Modes: LaplaceX / truth |
|---|---:|---:|---:|---:|
| single_uniform | 0.37 | 2.01 | 5.38× | 1 / 1 |
| two_uniform | 2.51 | 11.51 | 4.60× | 2 / 2 |
| three_log_pilot | 1.28 | 36.44 | 28.37× | 3 / 3 |
| two_log_pilot | 0.95 | 6.55 | 6.93× | 2 / 2 |
| three_log_holdout | 1.25 | 15.98 | 12.77× | 3 / 3 |
| close_rates_holdout | 7.99 | 50.77 | 6.35× | 2 / 2 |
| weak_fast | 2.41 | 18.31 | 7.61× | 3 / 3 |
| weak_middle | 2.11 | 13.33 | 6.32× | 3 / 3 |
| weak_slow | 2.09 | 40.00 | 19.13× | 3 / 3 |
| noisy_log_seed7 | 2.79 | 37.69 | 13.52× | 3 / 3 |
| noisy_log_seed19 | 2.58 | 36.42 | 14.13× | 3 / 3 |
| noisy_log_seed41 | 2.48 | 34.46 | 13.88× | 3 / 3 |

Times are medians of five local runs after warmup. SciPy comparison includes
both TRF and LM automatic-count fits with BIC, plus the original residual-stop
baseline. The table chooses the fastest configuration that recovers the count
and meets the reported reconstruction/rate accuracy criteria. Known-count
SciPy oracle results, every configuration, and individual timings are in the report.

### Where the accuracy advantage comes from

| Weak mode | LaplaceX modes / max abs error | SciPy residual-stop modes / max abs error | SciPy BIC / known-K accuracy |
|---|---:|---:|---|
| weak_fast | 3 / 4.44e-16 | 2 / 4.25e-04 | Also recovers 3 modes; error down to 6.66e-16 |
| weak_middle | 3 / 4.44e-16 | 2 / 3.79e-04 | Also recovers 3 modes; error down to 4.44e-16 |
| weak_slow | 3 / 5.20e-16 | 2 / 8.53e-04 | Also recovers 3 modes; error down to 4.44e-16 |

At 1000:1 amplitudes, the residual-stop SciPy baseline can stop at two modes.
Wavelet-guided selection recovers the weak third mode. Stronger SciPy BIC
selection or supplying the count also recovers it, with additional fitting
work. Successful fits generally reach the same float64 accuracy; this is
an automatic model-selection advantage, not universal precision superiority.

### Isolating the wavelet contribution

The pair-case follow-up caches noise thresholds at crossings, pads filter
boundaries, refines unresolved ridge proposals on both sides, and reuses QR
across damping trials. Count-selection rules and fit tolerances are retained.

| Previously slower case | Before: wavelet / grid ms | Now: wavelet / grid ms | Wavelet speedup vs previous | Now: grid / wavelet |
|---|---:|---:|---:|---:|
| two_uniform | 4.03 / 3.63 | 2.51 / 3.77 | 1.61× | 1.50× |
| close_rates_holdout | 12.87 / 8.68 | 7.99 / 8.23 | 1.61× | 1.03× |

A separate **21-run** randomized check of the final binary confirms the local pair timings.
Cells show median [25th, 75th percentile] milliseconds; every fit must pass the component-count and accuracy checks.

| Case | Default wavelets, ms | Rate grid, ms | Grid / wavelet |
|---|---:|---:|---:|
| two_uniform | 3.16 [2.95, 3.40] | 4.36 [4.10, 4.88] | 1.38× |
| close_rates_holdout | 7.14 [6.84, 7.50] | 8.33 [7.90, 8.75] | 1.17× |

The no-wavelet ablation uses rate-grid candidates with the **same new QR/LM
optimizer, noise tests and model selection**. It separates wavelet initialization
from the optimizer rewrite. Full before/after and ablation results are published.

With the same optimizer, default wavelet initialization is **19.25×** faster
on the widely separated three-mode case and **7.83×** faster on its holdout.
The two-uniform and close-rate ratios are **1.50×** and
**1.03×** respectively; above 1 means wavelets are faster.
These ratios compare the native rate-grid route
with the public wavelet route, including its Python routing overhead.

The exact single-mode shortcut is shared by both routes; speed gains on
that case belong to the shortcut rather than to wavelet analysis.

![Inverse speed and wavelet ablation](benchmarks/comparison/comparison.png)

The ordinary forward quadrature and Weeks inverse paths remain separate from
wavelet decomposition. Their speedups over mpmath are reported independently.
No finite benchmark establishes superiority on all signals: colored noise,
unobserved time scales and nearly coincident modes remain difficult.

[Full report, settings, historical measurements, and raw data](benchmarks/COMPARISON.md).

```bash
python -m pip install ".[benchmark]"
python benchmarks/validate_wavelets.py
python benchmarks/validate_wavelets.py --seed 20261107 --output benchmarks/comparison/validation_fresh.json
python benchmarks/validate_wavelets.py --seed 20261119 --output benchmarks/comparison/validation_holdout.json
python benchmarks/compare.py --section all --repeats 5 --points 64
python benchmarks/compare_extended.py --repeats 5 --points 64
python benchmarks/compare_pairs.py --repeats 21
python benchmarks/render_comparison.py --update-readme
```

<!-- comparison:end -->

## Quick Start

```python
import numpy as np
from laplacex import laplace, inverse, decompose, reconstruct, utils

# --- 1. Forward transform ---------------------------------------------------
# F(s) for f(t) = exp(-2t)  →  should equal 1/(s+2)
s_vals = np.array([1+0j, 2+0j, 3+1j])
F = laplace(lambda t: np.exp(-2*t), s_vals)

# --- 2. Inverse transform  ---------------------------
t = np.linspace(0.1, 5, 100)
f_rec = inverse(lambda s: 1/(s+2), t)          # recovers exp(-2t)

# --- 3. Decompose a noisy multi-exponential signal --------------------------
t   = np.linspace(0.01, 15, 600)
f   = utils.gen_exp_sum(t, [2.0, 0.5], [0.5, 5.0])
f_n = utils.add_noise(f, sigma=0.02, seed=0)

ec = decompose(t, f_n)       # → ExponentialComponents
print(ec)
# ExponentialComponents(2 components, significances=[0.832, 0.168])

for c in ec:
    print(f"  τ = {c.decay_time:.3f}   A = {c.amplitude:.3f}")

# --- 4. Reconstruct ---------------------------------------------------------
f_hat = reconstruct(ec, t)

# --- 5. Plot (requires matplotlib) -----------------------------------------
ec.plot(t_eval=t, t_original=t, f_original=f_n)
```

---

## API Reference

### Free functions

| Function | Description |
|---|---|
| `laplace(f, s, t=None)` | Forward Laplace transform F(s) |
| `inverse(F, t_eval, ...)` | Inverse Laplace transform f(t) |
| `decompose(t, f, cfg=None)` | Wavelet decomposition into `ExponentialComponents` |
| `reconstruct(comps, t_eval)` | Reconstruct f(t) from components |

### Classes

#### `LaplaceTransform`

Stateful engine — configure once, call many times. Caches decomposition
results so that `inverse(None, ..., method='auto')` uses the component
representation. An explicitly supplied transform callable takes precedence.

```python
from laplacex import LaplaceTransform

lt = LaplaceTransform(tol=1e-12, weeks_N=128)
F  = lt.laplace(lambda t: np.exp(-t), s=2+1j)
ec = lt.decompose(t, f)             # cached internally
f_rec = lt.inverse(None, t_eval, method='auto')
```

#### `ExponentialComponents`

Container for decomposition results with NumPy arrays, filtering, export,
and visualisation.

```python
ec.amplitudes    # ndarray, shape (n,)
ec.rates         # ndarray, decay rates λ_k
ec.decay_times   # ndarray, τ_k = 1/λ_k
ec.significances # ndarray, normalized absolute amplitudes

ec.filter(min_significance=0.1)
ec.top(3)
ec.reconstruct(t_eval)
ec.reconstruct_laplace(s_eval)
ec.to_dict()
ec.to_dataframe()   # requires pandas
ec.plot(t_eval, t_original, f_original)
```

### Configuration

```python
from laplacex import TransformConfig, DecomposeConfig

cfg = TransformConfig()
cfg.tol        = 1e-12     # absolute quadrature tolerance
cfg.rel_tol    = 1e-10
cfg.max_levels = 25
cfg.t_max      = 1e8

dcfg = DecomposeConfig()
dcfg.n_scales      = 16    # wavelet scales (default)
dcfg.n_log_points  = 256   # log-plane initialization grid (default)
dcfg.scale_min     = 0.05
dcfg.scale_max     = 0.5
dcfg.wavelet       = "mexh"
dcfg.snr_threshold = 5.0
dcfg.model_selection = "bic"  # or "aicc"
dcfg.use_aic       = True   # legacy name: enable information-criterion pruning
```

### Utilities (`laplacex.utils`)

```python
from laplacex import utils

t = np.linspace(0.01, 10, 500)
f = utils.gen_exp_sum(t, [1.0, 0.5], [1.0, 3.0])
f = utils.add_noise(f, sigma=0.05, seed=42)
s = utils.estimate_noise(f)
x, ft = utils.to_log_plane(t, f, n=512)
```

---

## Theoretical Background

Logarithmic coordinates `x = ln(t)` help analyse signals with widely separated
decay times. This does not remove the ill-conditioning of inverse Laplace
transforms or guarantee recovery of components from noisy data.

The decomposition algorithm:

1. Resamples `g(x) = f(exp(x))` with constant boundary extension.
2. Finds persistent Mexican-Hat zero crossings at extrema of the smoothed log derivative.
3. Tracks ridges from fine to coarse scales, retaining branches when they merge.
4. Fits all persistent components together on the original samples using an
   analytic Jacobian, Levenberg–Marquardt steps and scaled Householder QR.
5. Searches residual ridges for weak components and propagates sample noise
   through both interpolation and wavelet weights.
6. Selects component count using BIC by default, with uncertainty screening and
   refitted deletions; uses rate-grid proposals if the residual warrants a retry.

For `f(t) = A exp(-rate*t)`, the negative log derivative is
`A exp(z-exp(z))`, where `z=x+log(rate)`. Rates translate the same profile along
the log axis while preserving its height. This lets one wavelet scale range
initialize modes across several orders of magnitude. Log-grid estimates are
initial guesses; final rates and amplitudes are fitted to the original observations.

Only the Mexican Hat (`"mexh"`) wavelet is currently implemented.
Decomposition fits real decaying exponentials with signed amplitudes and recognizes
constant-only signals; it does not estimate oscillation frequencies. Explicitly supplied
oscillatory components can still be reconstructed in time and Laplace space.

## Numerical behaviour

Callable forward transforms integrate from zero to `cfg.t_max`. The tolerance
controls quadrature error, not the omitted tail: increase `t_max` when the
signal decays slowly. Failed adaptive convergence raises `RuntimeError`.

Sampled transforms integrate a natural cubic spline up to
`min(t[-1], cfg.t_max)`. If the first sample is after zero, the first spline
segment is extrapolated to zero. Samples imply zero signal beyond the final
time; include enough of the tail in your data. Times must be finite, increasing,
and nonnegative (strictly positive for decomposition).

The inverse returns real signals and uses the actual Weeks method. `N` controls
the number of Laguerre terms, `sigma` selects the contour to the right of all
singularities, and `b` selects the scale. Zero requests automatic selection
(`sigma = 4 / max(t_eval)`, with a search for `b`). Highly oscillatory signals
can require larger `N`; check convergence by increasing it. Automatic choices
are heuristics and do not provide an error guarantee. Use explicit parameters
for growing signals or difficult transforms. At `t=0` the returned value is
the finite expansion's right-limit estimate, when that limit exists.

Forward frequencies, inverse times, and component evaluation arrays preserve
their shapes, including scalars and NumPy slices. Sample data are one-dimensional.

## Development and verification

```bash
python -m pip install ".[dev]"
python -m pytest
python benchmarks/accuracy.py --output benchmarks/local.json
```

Building from source requires a C++17 compiler. The accuracy script uses known
analytic transform pairs; installing `mpmath` enables 60-digit reference
evaluation. [The verification report](benchmarks/ACCURACY.md) includes before
and after measurements and compatibility changes.

---

## License

MIT © 2026 Vozmishchev Artem
