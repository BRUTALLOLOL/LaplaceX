# LaplaceX 1.2.0

**Python library for exponential-mixture fitting and numerical Laplace transforms.**

LaplaceX combines adaptive Gauss–Kronrod forward quadrature, Weeks–Laguerre
inversion, and exponential fitting seeded by wavelet analysis in logarithmic
time. Accuracy depends on the signal, sampling window, and configuration;
see [evaluation and limitations](benchmarks/research/REPORT.md).

Version 1.2.0 is a recent research release. Tests establish a checked software
snapshot; the experimental example and single-machine timings do not establish
broad experimental validation or long-term adoption.

| Operation | Implementation |
|---|---|
| Callable forward transform | Adaptive embedded G7–K15 quadrature |
| Sampled forward transform | Natural cubic spline, integrated at every knot |
| Inverse transform | Weeks–Laguerre expansion with automatic scale search |
| Exponential decomposition | Persistent log-wavelet ridges, separable QR/LM fitting, BIC |
| Python interface | NumPy arrays and a compiled C++ core via pybind11 |

---


**Requirements:** Python ≥ 3.10, NumPy ≥ 1.24.
Optional extras: `matplotlib` (plotting), `pandas` (DataFrame export).

---

<!-- comparison:start -->
## Component recovery: version 1.2.0

Version 1.2.0 improves count recovery while retaining the original model-selection
criteria. Wavelets propose rates on a logarithmic time axis, followed by
separable QR/LM refinement and BIC/AICc selection on the original samples.
It estimates real decay rates and amplitudes, not oscillation frequencies.

| Evaluation | 1.1.0 wavelets | 1.2.0 wavelets | Variable projection + BIC |
|---|---:|---:|---:|
| Exposed development sweep: correct counts | 708/960 | 738/960 | 737/960 |
| Exposed development: count **and** accuracy | 573/960 | 613/960 | 613/960 |
| Independent sweep: correct counts | 700/960 | 730/960 | 728/960 |
| Independent sweep: count **and** accuracy | 562/960 | 596/960 | 596/960 |

The independent seed was fixed before evaluation. Both sweeps cover six rate
ratios (1.1–10), four Gaussian noise levels (0.0001–0.05), two sampling grids,
and twenty realizations per cell. The 480 original mixtures retain 480/480
correct counts. The primary wavelet route disables the Python rate-grid retry;
public-default results agree on the independent sweep.

The result is a tie with this variable-projection implementation on the joint
count/rate/reconstruction target. A two-case count advantage does not establish
statistically significant or universal superiority. Correct counts can still
have large rate errors: the revised independent-sweep 95th percentile is 41.11%
and the worst is 977.83%. Closely spaced noisy modes remain difficult.

### Speed and tradeoffs

Nine shuffled timing rounds on twelve existing fixtures retain all accuracy
checks. On the non-pilot multi-component fixtures, revised wavelets are
**1.64–13.66× faster than native rate-grid proposals with the same optimizer**.
The close-pair case improves from 5.89 to 2.11 ms. These familiar fixtures have
been used during development; the complete tables include the pilot cases.

Four fixtures are slower than the previous wavelet implementation, by roughly
8–59%, while retaining the advantage over the tested alternatives. The much
larger gap to the Python/SciPy variable-projection baseline includes native-code
and implementation effects; it is not attributed entirely to wavelets.
The single-mode shortcut is shared by both native routes.

A balanced timing check selects one fixed observation from every independent-sweep
cell: 48 cases, three rounds, all outcomes retained. Summing per-case medians,
revised wavelets take **0.745 s**, versus **16.939 s** for the Python baseline
(**22.73×** less time), with the same 37/48 counts and 31/48 joint successes.
The matched native grid takes **0.339 s**, so it is **2.20× faster than wavelets**
on this harder subset. This limits claims about a wavelet-specific speed advantage.

[Full comparison, thresholds, timing distributions and raw outcomes](benchmarks/COMPARISON.md).
[Original 1.0.0 research evaluation](benchmarks/research/REPORT.md), including
linear-coordinate controls, matrix pencil, information bounds and experimental
decays, remains a separate frozen record. The quadrature and Weeks-inversion
routines are unchanged by this detection revision.
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
ec.significances # ndarray, normalized absolute amplitudes; not p-values

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
