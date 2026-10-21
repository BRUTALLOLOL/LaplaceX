# Measured fluorescence data (CC BY 4.0)

Original data: Stephan Kupfer, Alexander Tombrink, Benedikt Callies, Daniel
Kowalczyk, and Benedikt Wiedemann, *Supplementary Data: 4CzIPN Photosensitizes
the Cobaloxime-catalyzed Light-Driven Hydrogen Evolution Reaction*, v3,
22 June 2026, https://doi.org/10.5281/zenodo.20792957.

License: Creative Commons Attribution 4.0 International,
https://creativecommons.org/licenses/by/4.0/.

The two CSVs are unchanged copies of the unquenched delayed and prompt
fluorescence files inside `Emission Decay.zip`, renamed for portable paths.
Original member names, units, and hashes are in `provenance.json`.
Each file includes an instrument response column; the fluorescence data are
in the separate `Decay` column. Reuse must credit the original authors.

Preprocessing (tail gating, background estimate, four-bin averaging) is described
in `../EXPERIMENT_PROTOCOL.md` and is applied only in memory, preserving originals.
This is a demonstration on experimental observations, without verified
component counts or rates. The native unweighted model is not an IRF-convolved
photon-likelihood analysis.
