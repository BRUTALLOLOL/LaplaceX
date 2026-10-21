## Component recovery: scope and evidence

Decomposition fits real decaying exponentials and selects a count up to the
component budget. It estimates amplitudes and rates; oscillation frequencies
are not estimated. Log-time wavelet proposals initialize a joint fit to the
original samples. Forward quadrature and Weeks inversion are separate routines.

The expanded evaluation keeps the 1.0.0 numerical core fixed. Its main route
disables the rate-grid retry in every section; public-default retry is a
separate sensitivity check. Every comparator's outcomes and failures are retained.

| Evaluation | Log wavelets | Variable projection + BIC | Matrix pencil + refinement |
|---|---:|---:|---:|
| Reused validation/regression mixtures: correct counts | 480/480 | 479/480 | 235/240 |
| New noise/separation sweep: correct counts | 708/960 | 737/960 | 323/480 |
| New sweep: correct count **and** rate/reconstruction target | 573/960 | 613/960 | 280/480 |

Matrix-pencil denominators include only originally uniform grids. The reused
sets include cases that informed development; the new sweep protocol was frozen
before evaluation. It uses two components, six rate ratios (1.1–10), four Gaussian
noise levels (0.0001–0.05), two sampling grids, and twenty realizations per cell.
Neither result establishes reliability for arbitrary decay data.

### What wavelet initialization improves

An external harness compiles log-time, linear-coordinate, and rate-grid routes
with the **same native optimizer, selection rules, compiler, and wrapper**.
Both coordinate detectors target peaks of the negative log-time derivative;
their point/scale budgets and relative filter widths match.

On the separated three-mode B fixture, the log route takes **1.23 ms**, versus
**8.88 ms** for linear coordinates and **12.29 ms** for rate-grid proposals
(**7.22×** and **9.99×** ratios). On the close pair, **linear coordinates are
faster: 6.80 ms versus 7.96 ms**. Times are medians of seven shuffled runs on
one Intel i7-7700HQ/Windows machine. All twelve development-exposed fixtures,
timings, and accuracy checks are reported; the single-mode shortcut is shared.

These controls support a reduction in proposal-search cost for some multiscale
signals. The earlier SciPy timing ratios include Python callback and solver
differences and are not used to isolate the wavelet contribution.

### Count recovery and parameter uncertainty

Correct counts do not ensure correct individual rates. The new sweep's
95th-percentile rate error among correctly counted log-wavelet fits is **31.85%**,
with extreme failures retained in the raw records. The joint target above
requires maximum rate error ≤10% and clean-signal RMSE ≤noise standard deviation.
The old 21.07% worst noisy case is also reproduced by a true-count variable-
projection fit; a local information calculation indicates weak identifiability.
Cramér–Rao diagnostics assume known count and independent Gaussian noise and do
not certify automatic model selection.

The clean 1000:1 weak-mode examples distinguish a premature residual-stop rule
from stronger count selection: SciPy BIC and known-count fits also recover the
weak modes. Successful clean fits generally reach comparable float64 accuracy.

Two openly licensed measured fluorescence tails add an application example.
Their true component counts are unknown, and an additional gate delay changes
the delayed-fluorescence count. Instrument response, background, and photon-noise
models matter before assigning physical lifetimes.

[Expanded report and reproduction instructions](benchmarks/research/REPORT.md) ·
[Frozen protocol](benchmarks/research/PROTOCOL.md) ·
[Release comparisons and historical measurements](benchmarks/COMPARISON.md) ·
[Historical numerical audit](benchmarks/ACCURACY.md)
