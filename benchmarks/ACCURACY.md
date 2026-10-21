# Numerical verification — 20 October 2026

This report preserves the initial numerical fixes, before the subsequent
log-plane wavelet rewrite. For the current component-recovery results,
160-mixture validation, and SciPy comparisons, see [COMPARISON.md](COMPARISON.md).

LaplaceX was compiled before and after the changes and compared against analytic
references using the same inputs. This report describes these cases; it does not
claim uniform accuracy for arbitrary inverse transforms or component mixtures.

## Method

- Baseline: committed C++ and Python sources at `6ba2877`, compiled separately
  from an untouched source snapshot into `.build/baseline/laplacex`.
- Revised version: working-tree sources compiled into `src/laplacex`.
- Both comparison builds: Zig 0.16.0 / Clang 21.1.0, C++17, `-O2`,
  `-target x86_64-windows-gnu`; Python 3.14.7; NumPy 2.5.3; Windows x64.
- The inverse closed-form references were evaluated with mpmath at 60 decimal
  digits and rounded to float64 for comparison. Other analytic references use
  float64 arithmetic. Roundoff-level differences are not meaningful speed or
  accuracy guarantees.
- Inverse grids: 100 points from 0.1 to 5. Decomposition grids: 600 points
  from 0.01 to 15. Exact definitions are in [accuracy.py](accuracy.py).
- Raw measurements: [before.json](before.json), [after.json](after.json).

## Maximum absolute errors

| Calculation | Before | After |
|---|---:|---:|
| Forward exponential, four real/complex frequencies | 1.661e-15 | 7.083e-16 |
| Forward fast decay (rate and frequency `1e12`) | 9.990e-16 | 1.010e-28 |
| Forward `1/sqrt(t)` at `s=2` | 7.105e-08 | 2.220e-16 |
| Sampled exponential, 2,000 samples | 4.755e-07 | 2.236e-08 |
| Inverse `1/(s+2)`, default N=64 | 1.188e-05 | 1.609e-14 |
| Inverse `1/s`, default N=64 | 1.240e-05 | 3.397e-14 |
| Inverse `1/s^2`, default N=64 | 6.698e-06 | 1.066e-14 |
| Inverse damped cosine, default N=64 | 1.192e-01 | 3.288e-12 |
| Inverse `10/(s^2+100)`, N=64 | 1.035e+00 | 6.159e-02 |
| Inverse `10/(s^2+100)`, N=256 | 1.035e+00 | 3.538e-09 |
| Analytic oscillatory component image, nonzero phase | 7.793e-01 | 2.776e-17 |
| Single-exponential reconstruction | 5.593e-08 | 7.261e-14 |
| Two-exponential reconstruction | 1.248e-03 | 3.359e-12 |

The fast-decay integral is only `5e-13`; its relative error fell from
`1.998e-3` (about 0.2%) to `2.019e-16`. Ordinary exponential forward transforms
were already accurate near float64 roundoff; the main gain there is cost.

## Performance and tradeoffs

The four-frequency forward check decreased Python callback evaluations from **2,426,805**
to **1,695**, a **1432× reduction**. This is a deterministic count,
not a timing estimate. Correct embedded Gauss weights allow adaptive integration
of smooth intervals to stop instead of unnecessarily bisecting them.

Wavelet values are precomputed once per scale. Rate fitting reuses unchanged
basis columns and stores each column contiguously. Householder QR replaces
normal equations, avoiding the extra condition-number squaring. Rate refinement
now has convergence checks instead of stopping after three passes.

Raw wall times are single-run observations and can vary with machine load.
In this run, single-exponential fitting decreased from about 0.36 s to 0.025 s;
two-exponential fitting decreased from about 0.44 s to 0.31 s while fitting more
accurately. Weeks inversion costs more than the previous 18-term Stehfest
calculation for these grids. Its parameter search and additional contour samples
are the measured accuracy/performance tradeoff; the old implementation ignored
`sigma`, `b`, and orders above 18.

## Bugs corrected

- Embedded G7 weights were attached to the wrong K15 nodes.
- Callable quadrature omitted `[0, 1e-15]` and could exceed `t_max`.
- Integrable diffusion singularities now use a squared-coordinate endpoint
  substitution; nonconvergence and nonfinite callbacks raise explicit errors.
- The function labeled Weeks inversion actually implemented Stehfest. It now
  computes Weeks Fourier coefficients and scaled Laguerre recurrences, sharing
  coefficients across evaluation times.
- Oscillatory component transforms omitted the conjugate pole and phase.
- Array pointer reads ignored NumPy strides. Inputs now request contiguous
  storage before native conversion; evaluation shapes are preserved.
- Sample lengths, ordering, finite values, configurations, inverse orders, and
  inverse evaluation times are checked.
- A stale component cache no longer overrides an explicitly supplied image.
- AIC uses its sample-count factor and handles exact fits correctly.
- Release tests no longer exclude the test module or suppress failures.

## Compatibility and remaining limits

- Sampled forward transforms now integrate to `min(t[-1], cfg.t_max)`, using
  zero continuation beyond the final sample. The first spline segment is
  extrapolated to zero when the first observation is positive. This replaces
  the old constant extension at both ends; supply adequate endpoint samples.
- `t_max` truncation error is separate from quadrature tolerance. Slowly decaying
  signals need a sufficiently long window.
- The inverse assumes a real time-domain signal. Automatic `sigma` and `b` are
  heuristics, not certified error bounds. Set the contour to the right of all
  singularities for growing signals. Widely separated time scales, singular
  right limits, discontinuities, and strong oscillations need convergence checks.
- The sine example still has approximately 0.062 absolute error at N=64. N=256
  gives `3.54e-9`; increasing N now changes the actual expansion order.
- Exponential decomposition remains ill-conditioned for closely spaced rates
  and noisy mixtures. A separate clean probe with rates 0.5 and 0.7 and amplitudes
  2 and 0.5 over `[0.01,20]` recovered approximately 0.50161 and 0.71114, with
  normalized reconstruction error `3.36e-5`. Accurate reconstruction does not
  guarantee exact component identification.
- Only `mexh` is implemented; other wavelet names now raise a clear error.
- Native CPU tuning is opt-in (`-DLAPLACEX_NATIVE=ON`), so default wheels do
  not require the build machine's instruction set.

## Verification and reproduction

**99 tests passed** against both the local compiled extension and the installed
Windows CPython 3.14 wheel. The wheel was built through the normal
`scikit-build-core` / CMake / Ninja path, installed into `.venv`, and tested with
`python -I -m pytest` to avoid importing the source package accidentally.
The installed wheel reproduced the same numerical errors and callback counts.
Linux, macOS, and other Python versions were not executed locally; automated
Linux/Windows tests and release wheel tests are configured for those builds.

```bash
python -m pip install ".[dev]" mpmath
python -m pytest
python benchmarks/accuracy.py --output benchmarks/local.json
```

For a before/after comparison, build the baseline commit in a separate checkout
with the same compiler and flags, then run the same benchmark script against
each installed version. Preserve both JSON reports. A C++17 compiler is required.

## Numerical references

The Weeks formulation follows the Laguerre expansion described by
[Abate, Choudhury and Whitt (1996)](https://www.columbia.edu/~ww2040/Laguerre1996.pdf)
and the [University of Arizona Weeks method notes](https://archive.math.arizona.edu/brio/WEEKS_METHOD_PAGE/pkanoWeeksMethod.htex).
Parameter selection remains a separate numerical problem; see
[Weideman (1999)](https://epubs.siam.org/doi/10.1137/S1064827596312432).
