# Expanded research evaluation

The numerical core is frozen at release 1.0.0. These local research records do not imply a new package release.
The [protocol](PROTOCOL.md) was fixed before inspecting outcomes. No core algorithm was tuned to the new sweep.

## Count recovery

| Method | Legacy correct / eligible | New sweep correct / eligible | Sweep count + accuracy pass |
|---|---:|---:|---:|
| Log wavelets | 480/480 | 708/960 | 573/960 |
| Linear-coordinate wavelets | 477/480 | 708/960 | 580/960 |
| Native rate grid | 480/480 | 737/960 | 568/960 |
| Variable projection + BIC | 479/480 | 737/960 | 613/960 |
| Pencil + variable projection | 235/240 | 323/480 | 280/480 |

Legacy cases are reused validation/regression data. New sweep cells have 20 noise realizations each.
Rate errors are conditional on correct counts; the last column counts all eligible observations, including incorrect counts as failures.
Pencil methods are evaluated only on originally uniform samples. Known-K oracles are not scored as automatic count recovery.

## Matched native timing

| Fixture | Log ms | Linear ms | Grid ms | VarPro ms | Grid/log | Linear/log |
|---|---:|---:|---:|---:|---:|---:|
| single_uniform | 0.155 | 0.192 | 0.163 | 8.169 | 1.05 | 1.24 |
| two_uniform | 2.895 | 2.968 | 4.398 | 26.218 | 1.52 | 1.03 |
| three_log_pilot | 1.364 | 21.004 | 30.387 | 274.506 | 22.27 | 15.40 |
| two_log_pilot | 0.833 | 2.411 | 4.298 | 22.048 | 5.16 | 2.90 |
| three_log_holdout | 1.229 | 8.880 | 12.285 | 64.593 | 9.99 | 7.22 |
| close_rates_holdout | 7.963 | 6.801 | 9.783 | 36.819 | 1.23 | 0.85 |
| weak_fast | 2.111 | 194.037 | 11.200 | 81.569 | 5.31 | 91.91 |
| weak_middle | 1.952 | 7.457 | 9.770 | 48.922 | 5.01 | 3.82 |
| weak_slow | 1.994 | 7.414 | 10.109 | 64.836 | 5.07 | 3.72 |
| noisy_log_seed7 | 2.665 | 5.765 | 13.692 | 236.877 | 5.14 | 2.16 |
| noisy_log_seed19 | 2.191 | 4.521 | 11.982 | 146.059 | 5.47 | 2.06 |
| noisy_log_seed41 | 2.146 | 5.824 | 13.483 | 488.796 | 6.28 | 2.71 |

Seven shuffled repetitions after warmup; one BLAS thread. Every repetition and its accuracy is retained in timing.json.
A ratio above one favors the log route. Fast failed fits are not accuracy-preserving speedups.
The single-mode shortcut is shared. Existing pilot/holdout fixture names do not indicate untouched data.
The coordinate control changes both finite-grid allocation and the geometry of smoothing while targeting the same -t*f'(t) peaks.
It is not a claim that all linear-time wavelet methods must fail or that complexity is logarithmic.

## Noise, separation, and local information

Rates are 0.5 and 0.5*r, amplitudes 1 and 0.5. Ratios 1.1, 1.25, 1.4, 2, 3, 10; noise SD 0.0001, 0.001, 0.01, 0.05.
Both sampling grids use 600 observations on [0.01,20]. Seed 2026102101. All 960 outcomes are retained.
CR bounds assume known count, known independent Gaussian noise, and local unbiased estimation. Amplitudes are nuisance parameters.
Known-K empirical log-rate RMSE and bias use all 20 realizations per cell, without conditioning on automatic count success.
A single error, a biased estimator, and a selected model are not governed directly by the displayed local standard deviations.

The legacy worst-rate-error case is 20261107-uniform-55. Native and known-K variable projection reproduce the same poorly determined fit.
## Experimental application

Two unquenched TCSPC tails from [Kupfer et al., Zenodo v3](https://doi.org/10.5281/zenodo.20792957), CC BY 4.0.
Original files and hashes are in [data](data/README.md). [Preprocessing rules](EXPERIMENT_PROTOCOL.md) were fixed before fitting.
Component counts and rates are not known. Gate sensitivity and unequal residuals are reported, not treated as recovery accuracy.
Photon noise is heteroscedastic and the IRF is gated rather than convolved. This is a descriptive application outside the Gaussian validation model.
The rounded experimental time columns fail the strict uniform-grid eligibility check; pencil fits are not included for those curves.

## Configuration sensitivity

- Public default (retry enabled): 480/480 legacy counts. This is a secondary configuration, not a replacement for the main row.
- Linear coordinates, 1024 points: 376/480 legacy counts. This is a secondary configuration, not a replacement for the main row.

The 1024-point linear detector performs worse here; increasing resolution is not a recovery guarantee.
## Reproducibility

CPU: Intel Core i7-7700HQ, nominal 2.80 GHz; Balanced power plan. No fixed-clock or thermal-control claim.
Compiler and source hashes, Python/NumPy/SciPy versions, exceptions, warnings and raw fits are in results/.
The external harness compiles a copy of cpp/src/critical.cpp and a checked coordinate transform; it does not replace the installed extension.
Independent analytic-Jacobian, Fisher-information, baseline, and native-equivalence checks are included.

```powershell
.venv\Scripts\python benchmarks\research\build_native.py --compiler .venv\Lib\site-packages\ziglang\zig.exe
.venv\Scripts\python benchmarks\research\check_native.py
.venv\Scripts\python benchmarks\research\run_study.py --section all
.venv\Scripts\python benchmarks\research\experimental.py
.venv\Scripts\python benchmarks\research\summarize.py
```
