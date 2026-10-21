# Component-detection revision protocol

The 1.1.0 core and the research baselines are the before-change controls.
The original 960-signal sweep is development data for this revision. Never
label it untouched validation or overwrite benchmarks/research/results.

Primary endpoints: correct component count and joint count/rate/reconstruction
success, using the previous thresholds. Report conditional rate errors as well;
count success alone can conceal poor rates. Compare the actual automatic-count
variable-projection implementation, with its original two starts and BIC.

Validation after selecting the implementation: repeat the 960-case protocol
with seed 2026102207 (not used for development), plus the existing 480 mixtures.
Keep all failures. Evaluate the native wavelet route with retry disabled and
the public default separately. Also report noise/constant and signed-mixture
regressions. No algorithm tuning after inspecting this validation seed; if
further tuning is needed, mark that run development and use another fresh seed.

Speed: the same twelve archived timing fixtures, at least seven shuffled runs
after warmup, one BLAS thread, before and after in the same process. Include
matched native grid and automatic-count variable projection, all failures and
rate/count accuracy. Report per-fixture regressions as well as aggregate gains.

Do not modify or upload manuscript files. Update benchmarks/COMPARISON.md with
the new measurements and limits before requesting approval to commit or push.

## Candidate revision before the second validation run

The first candidate was rejected after the signed-mixture supplement exposed
missed clean four-mode signals. Its partially completed seed-2026102207 run
(660 observations) is development data, not independent validation. Its raw
records and source snapshot are retained locally under .build/detection-rejected-1.
The screen had treated model mismatch as measurement noise. The second candidate
caps its screening variance by the measured-noise estimate as well as fitted
residual variance. A second signed-mixture failure came from inadmissible
coalescing-rate projected fits with unbounded cancelling amplitudes. Only when
the projected fit violates the pre-existing admissibility limits, retry the
original joint optimizer from the original proposals. Accept an admissible
joint result rather than letting a forbidden projected minimum block growth.
The selected implementation restores both clean four-component observations.

The final independent sweep now uses seed 2026102217, with the same 960-case
layout and endpoints. This addendum precedes that run. No tuning after inspecting
that seed is permitted without again declaring it development data and reseeding.
