# Signed-mixture and null-signal checks

Fixed before running: RNG seed 2026102209; grids uniform [0.01,20] and
geometric [0.001,1000], 600 samples. Counts 1–4, ten observations per count,
noise standard deviations 0, 0.001, 0.01: 240 signed mixtures. Draw absolute
amplitudes log-uniformly from 0.2–2 and independent signs. Draw rates with the
same bounds/separation as the original 480-mixture protocol. Component truth
is not changed after cancellation or noise is observed. Apply the same joint
accuracy targets; retain all failures.

Also test twenty independent noise observations per grid with standard
deviation 0.001, and the same observations plus a unit constant: 80 null or
constant checks. A constant is one zero-rate component. Use the before/after
native routes. This supplement checks regressions, not a claim of complete
coverage of all signed or non-white-noise signals.
