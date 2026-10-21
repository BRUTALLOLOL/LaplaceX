"""Fixed independent validation protocol; no SciPy timing claims use this sample.

Final validation seed 20261029; 160 mixtures, 600 points each, 1--4 components; half noisy.
Log grids: rates 0.02--60, neighboring rate ratio >=3, t=0.001--1000.
Uniform grids: rates 0.2--6, neighboring rate ratio >=2, t=0.01--20.
Amplitudes 0.2--2, log-uniform; white Gaussian noise sigma=0 or 0.001.
The fallback is disabled to test component recovery by the wavelet route itself.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np

import laplacex as lx
from laplacex import _laplacex as native


def run(seed=20261029):
    rng = np.random.default_rng(seed)
    rows = []
    for grid in ("log", "uniform"):
        for index in range(80):
            count = 1+index%4
            low, high, separation = (.02, 60., 3.) if grid == "log" else (.2, 6., 2.)
            while True:
                rates = np.sort(np.exp(rng.uniform(np.log(low), np.log(high), count)))
                if count == 1 or min(rates[1:]/rates[:-1]) >= separation:
                    break
            amplitudes = np.exp(rng.uniform(np.log(.2), np.log(2), count))
            t = np.geomspace(.001, 1000, 600) if grid == "log" else np.linspace(.01, 20, 600)
            clean = np.exp(-t[:, None]*rates)@amplitudes
            sigma = 0. if index < 40 else .001
            observed = clean+rng.normal(0, sigma, len(t))
            start = perf_counter()
            result = lx.decompose(t, observed, auto_deflate=False)
            elapsed = perf_counter()-start
            recovered = np.sort(result.rates)
            count_ok = len(recovered) == count
            rate_error = float(max(abs(recovered/rates-1))) if count_ok else None
            error = result.reconstruct(t)-clean
            rows.append({
                "grid": grid, "index": index, "rates": rates.tolist(),
                "amplitudes": amplitudes.tolist(), "noise_sigma": sigma,
                "recovered_rates": recovered.tolist(), "correct_count": count_ok,
                "max_rate_relative_error": rate_error,
                "max_abs_error": float(max(abs(error))),
                "clean_rmse": float(np.sqrt(np.mean(error**2))), "seconds": elapsed,
            })
    return {
        "seed": seed, "protocol": __doc__, "module": lx.__file__,
        "extension_sha256": hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "correct_counts": sum(row["correct_count"] for row in rows),
        "cases": len(rows), "results": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent/"comparison/validation.json")
    args = parser.parse_args()
    result = run(args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(f"Correct component counts: {result['correct_counts']}/{result['cases']}")
    for row in result["results"]:
        if not row["correct_count"]:
            print(row)
