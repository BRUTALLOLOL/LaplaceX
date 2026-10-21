# Additional balanced timing check

Before timing, select repeat 0 from every cell in the completed independent
seed-2026102217 sweep: 2 grids × 6 rate ratios × 4 noise levels = 48 cases.
Do not select by correctness, residual, or elapsed time. Repeat the before/after,
matched native-grid and unchanged variable-projection fits in three shuffled
round-robin rounds after warmup, with one BLAS thread and no competing benchmark
jobs. Retain incorrect counts and parameter errors in every timed record; no
fastest-passing-comparator selection. This is supplementary latency coverage,
not another independent accuracy sample or a new tuning dataset.
