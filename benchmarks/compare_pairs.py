"""Longer same-run wavelet/rate-grid check for the two formerly slower cases."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from compare import HERE, Output, decomposition_cases, lx, measure, metadata, native


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=21)
    parser.add_argument("--output", type=Path, default=HERE/"comparison/pairs.json")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    args.points = 64
    results = {}
    for name, t, clean, observed, rates, amplitudes, sigma, seed in decomposition_cases():
        if name not in ("two_uniform", "close_rates_holdout"):
            continue

        def components(output):
            if not isinstance(output, lx.ExponentialComponents):
                output = lx.ExponentialComponents(output)
            return Output(output.reconstruct(t), output.rates)

        cfg = lx.DecomposeConfig()
        methods = {
            "LaplaceX default": lambda: components(lx.decompose(t, observed)),
            "LaplaceX wavelets only": lambda: components(lx.decompose(t, observed, auto_deflate=False)),
            "LaplaceX no wavelets": lambda: components(native.decompose_deflate(t, observed, cfg, 5, 1e-6)),
        }
        results[name] = {
            "rates": rates.tolist(), "amplitudes": amplitudes, "noise_sigma": sigma,
            "methods": measure(methods, clean, args.repeats, observed, rates),
        }
    meta = metadata(args)
    meta["pair_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"metadata": meta, "results": results}, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print("Saved 21-run pair check" if args.repeats == 21 else f"Saved {args.repeats}-run pair check")
