# Response to the scientific and repository review

The revision retains the 1.0.0 numerical core and adds external evaluation rather
than tuning the algorithm to the new data. The companion evaluation and documentation are included in version 1.1.0.

| Review point | Revision and evidence | Remaining limit |
|---|---|---|
| C++/Python timing confound | Same compiled wrapper, compiler, QR/LM optimizer and selection code for log, linear-coordinate, and grid proposals; twelve fixtures timed seven times | One machine; specialist Python timing still includes callback overhead |
| No coordinate ablation; ambiguous title | Title concerns fitting mixtures with wavelet proposals. Both coordinate controls target peaks of -t*f'(t), with matched budgets | Tests allocation and smoothing geometry together; better linear-time detectors may exist |
| 480 cases lack comparators | Variable projection on all 480; pencil + refinement on the 240 eligible uniform cases; errors and counts retained | These remain reused validation/regression data |
| Easy noise/separation regime | Frozen sweep of 960 additional signals, six rate ratios, four noise levels, twenty realizations per cell | Two positive components and independent Gaussian noise only |
| No information bound | Full amplitude/log-rate information, nuisance amplitudes, rank checks, known-count Monte Carlo RMSE and bias | Local unbiased known-count assumptions; no selected-model guarantee |
| Thin related work | 18 references including Alfrey, Schwarzl--Staverman, Istratov--Vyvenko, Prony, HSVD, matrix pencil, NNLS and CONTIN | Integration/software contribution; no priority claim |
| Mixed configurations | Primary route: default detector, maximum five modes, retry disabled throughout. Public retry and 1024 linear points are separate | Public-default sensitivity uses reused data |
| No measured data | Two CC BY fluorescence tails, original files and hashes, predeclared preprocessing and gate sensitivity | Unknown component truth; no IRF-convolved photon-likelihood model |
| Statistical-sounding attribute | Local preferred `weights` alias, compatibility retained, two added tests | Not yet released to PyPI |
| README framing and development log | Condensed current evidence, distinct accuracy/count/speed claims, historical audit labelled | Live GitHub description awaits permission to change |
| Recent code and earlier serious bugs | Maturity stated, audit preserved, published CI and publishing verified | No long-term adoption or independent validation claimed |

## Findings that changed the paper

The stronger evaluation does not establish general supremacy. On a separated
three-mode fixture log proposals take 1.23 ms versus 8.88 ms for the matched
linear detector. The linear control is faster on the close pair. On the new
sweep both coordinate detectors recover 708/960 counts; variable projection
recovers 737/960. Joint count/accuracy success is 573/960 versus 613/960.

The earlier worst noisy rate error is reproduced by a true-count fit, with
a local log-rate CR standard deviation of 0.197 for that rate. This supports
a weak-identifiability interpretation without proving optimality. The harder
sweep also contains extreme rate errors despite correct counts. Fluorescence
gate sensitivity is reported as a limitation.

These changes strengthen a scoped software paper. They do not create the broad
reliability, independent adoption, or mathematical novelty missing previously.
