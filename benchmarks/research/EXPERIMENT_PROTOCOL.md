# Experimental example: preprocessing fixed before fitting

Source: Kupfer et al., Zenodo 20792957 (version v3, 22 June 2026),
DOI 10.5281/zenodo.20792957, CC BY 4.0. Archive: Emission Decay.zip,
MD5 65718629f717491b3c043f1f87e22b16.

Use the unquenched files ending S1-DF.csv and S1-PF.csv. These contain measured
time, instrument response (IRF), and fluorescence decay in separate columns.
Read the named Decay column, not the IRF column. Preserve original time units
(microseconds for DF, nanoseconds for PF). No component count or rate truth is
available for these curves.

For each curve, start the tail after the later of (i) five bins beyond the
measured decay's maximum and (ii) the last IRF bin exceeding 1% of its maximum.
Estimate a constant background as the median of the final 10% of decay bins.
Average non-overlapping groups of four consecutive tail bins, discarding only
an incomplete final group. Subtract the estimated background without clipping
negative values. Retain the original positive time coordinates. Do not tune
the gate, background rule, or bin size from fit outcomes.

Fit all the retained bins with the same native log/linear/grid and automatic
variable-projection and matrix-pencil algorithms. Report fitted counts, rates,
observed RMSE and residuals. Refit with an additional five-bin gate delay as a
predeclared gate sensitivity. No count or parameter 'accuracy' is asserted.
Unweighted fits are descriptive: photon-counting noise is heteroscedastic, the
instrument response is excluded by tail gating rather than deconvolved, and
the remaining background and model inadequacy can affect fitted rates.
