# Frozen evaluation protocol, 21 October 2026

This protocol is written before inspecting results. The released decomposition
source at ff164d309a89b72eaafb2c58eb83430399ba9b24 is held fixed. Failures are
reported; no detector or model-selection tuning follows examination of this run.
This file's SHA-256 is recorded in the output.

## Configurations

The main LaplaceX accuracy configuration is the native log-wavelet route, maximum
five components, default detector settings, no rate-grid retry. Public default
results are a separate sensitivity analysis. Both are evaluated on the same data.
The external native harness compiles an unchanged copy of critical.cpp for the
log and rate-grid routes. A second copy changes only the detector's sampling
coordinate, derivative geometry, and conversion from peak position to rate.
It shares the optimizer, selection, noise handling, and component budget.

## Specialist comparators

- Variable projection: exact QR-based projected-residual Jacobian; signed linear
  amplitudes; positive bounded rates; two deterministic starts for each count
  1--5; BIC selection. Both automatic-count and true-count results are retained.
- Matrix pencil: truncated-SVD Hankel pencil, followed by variable-projection
  refinement; counts 1--5 selected by BIC. Only originally uniform samples are
  eligible. Logarithmic samples are not interpolated to make the method applicable.
- Native rate-grid route: compiled with the same compiler and wrapper as the
  wavelet routes. This remains the implementation-matched non-wavelet control.

## Existing 480 mixtures

Regenerate observations exactly from the three archived seed protocols and verify
their true rates/amplitudes against every stored row. Evaluate all eligible
comparators, preserving every count, parameter error, reconstruction error,
elapsed time, and exception. Report errors separately from count recovery.
These are reused data, not a new independent validation set.

## Untouched noise and separation sweep

Two positive modes with rates 0.5 and 0.5*r, amplitudes 1 and 0.5, 600 samples.
Rate ratios r = 1.1, 1.25, 1.4, 2, 3, 10. White Gaussian standard deviations
0.0001, 0.001, 0.01, 0.05. Uniform or geometric grids on [0.01, 20]. Twenty
independent noise realizations per cell, RNG seed 2026102101: 960 observations.
Evaluate the same automatic-count methods plus true-count variable projection.
No fit is discarded because it fails count recovery or exceeds an error cutoff.

For each cell compute the Gaussian known-count Cramer--Rao covariance from the
full amplitude/log-rate Jacobian at the true parameters. Report rank/conditioning;
do not substitute a pseudoinverse's zero null-space variance. Compare unconditioned
true-count empirical log-rate RMSE and bias against local CR standard deviations.
These bounds do not govern biased automatic model selection or a single error.

## Coordinate ablation and timing

Use the twelve existing fixtures plus the 960 new sweep observations. Both
detectors target extrema of -t*f'(t). The log detector smooths f(exp(x)); the
linear detector smooths f(t) and uses the chain-rule factor t. Both use 256
points, 16 scales, the same coordinate span and relative filter widths. Linear
coordinates are an affine mapping of t onto the same span as log(t). This tests
fixed-budget allocation and smoothing geometry together, not an asymptotic claim.
Record initial proposals as well as final counts. Also evaluate 1024 linear
points as a separate budget sensitivity, never as the main matched-budget row.

Time native routes and variable projection over seven randomized round-robin
repetitions after a warmup on all twelve fixtures. Report the two pilot fixtures
separately; do not choose one passing comparator per fixture for the new headline.
Use one BLAS thread. Report all outcomes even when the fastest method fails.

## Experimental example

Use a publicly archived measured decay with identifiable provenance and license.
Record the original file hash and all tail-selection/background rules before
fitting. Compare reconstructions and residuals. Without independent component
truth, report fitted rates as descriptive parameters, not verified molecular
components. Photon-counting noise and instrument-response effects must be stated.
