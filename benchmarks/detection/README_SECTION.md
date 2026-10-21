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
