"""Reproducible comparisons and wavelet ablations; no library code is modified.

python benchmarks/compare.py --section all --repeats 5 --points 64
Requires numpy, scipy, mpmath, matplotlib, and the compiled LaplaceX package.
The fine-wavelet preset was explored on three pilot signals; the report also
includes default settings, failures, and additional signals rather than hiding them.
"""

from __future__ import annotations

import os

# Apply before importing NumPy/SciPy so small fits do not oversubscribe BLAS.
for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import argparse
import cmath
import hashlib
import json
import math
import platform
import random
import statistics
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter

import mpmath as mp
import numpy as np
import scipy
from scipy import integrate, optimize, signal

import laplacex as lx
from laplacex import _laplacex as native

HERE = Path(__file__).resolve().parent


@dataclass
class Output:
    values: np.ndarray
    rates: np.ndarray | None = None
    metadata: dict = field(default_factory=dict)


def summarize(output, truth, observed=None, true_rates=None):
    values = np.asarray(output.values)
    if values.shape != truth.shape or not np.all(np.isfinite(values)):
        raise ValueError("Output shape mismatch or nonfinite output")
    error = abs(values-truth)
    result = {
        "max_abs_error": float(np.max(error)),
        "normalized_max_error": float(np.max(error)/np.max(abs(truth))),
        "rmse": float(np.sqrt(np.mean(error**2))),
        "passes_transform_target": bool(np.all(error <= 1e-10+1e-8*abs(truth))),
        **output.metadata,
    }
    if observed is not None:
        result["training_relative_l2"] = float(np.linalg.norm(values-observed)
                                               /np.linalg.norm(observed))
    if true_rates is not None:
        recovered = np.sort(output.rates)
        expected = np.sort(true_rates)
        result["rates"] = recovered.tolist()
        result["component_count"] = len(recovered)
        result["correct_component_count"] = len(recovered) == len(expected)
        result["max_rate_relative_error"] = (
            float(np.max(abs(recovered/expected-1))) if len(recovered) == len(expected)
            else None
        )
    return result


def measure(methods, truth, repeats, observed=None, true_rates=None):
    """One untimed warmup, then randomized round-robin timed runs; include setup."""
    records, active = {}, {}
    for name, function in methods.items():
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                output = function()
            records[name] = {
                "status": "ok",
                "warnings": sorted({str(w.message) for w in caught}),
                **summarize(output, truth, observed, true_rates),
                "timings_seconds": [],
            }
            active[name] = function
        except Exception as error:
            records[name] = {"status": "error", "error": f"{type(error).__name__}: {error}"}
    rng = random.Random(20261020)
    for _ in range(repeats):
        order = list(active)
        rng.shuffle(order)
        for name in order:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                start = perf_counter()
                output = active[name]()
                elapsed = perf_counter()-start
            records[name]["timings_seconds"].append(elapsed)
            records[name].update(summarize(output, truth, observed, true_rates))
    for item in records.values():
        if item["status"] == "ok":
            samples = item["timings_seconds"]
            item["median_seconds"] = statistics.median(samples)
            item["min_seconds"] = min(samples)
            item["max_seconds"] = max(samples)
    return records


def reference(function, t):
    with mp.workdps(60):
        return np.array([float(function(mp.mpf(float(x)))) for x in t])


def forward_suite(repeats, points):
    results = {}
    for name, function, frequencies, exact in [
        ("exp_four", lambda t: math.exp(-2*t), np.array([.1, 1, 3+1j, 10+20j]),
         lambda s: 1/(s+2)),
        ("exp_batch", lambda t: math.exp(-2*t), np.geomspace(.1, 100, points).astype(complex),
         lambda s: 1/(s+2)),
        ("diffusion_four", lambda t: 1/math.sqrt(t), np.array([.1, 1, 3+1j, 10+20j]),
         lambda s: np.sqrt(np.pi/s)),
    ]:
        truth = exact(frequencies)
        cfg = lx.TransformConfig()
        cfg.tol, cfg.rel_tol = 1e-10, 1e-8

        def quad_each():
            values = []
            for s in frequencies:
                if s.imag == 0:
                    value = integrate.quad(lambda t: function(t)*math.exp(-s.real*t),
                                           0, np.inf, epsabs=cfg.tol, epsrel=cfg.rel_tol,
                                           limit=300)[0]
                else:
                    value = integrate.quad(lambda t: function(t)*cmath.exp(-s*t),
                                           0, np.inf, epsabs=cfg.tol, epsrel=cfg.rel_tol,
                                           limit=300, complex_func=True)[0]
                values.append(value)
            return Output(np.asarray(values))

        methods = {
            "LaplaceX": lambda: Output(lx.laplace(function, frequencies, cfg=cfg)),
            "SciPy quad": quad_each,
            "SciPy quad_vec": lambda: Output(integrate.quad_vec(
                lambda t: function(t)*np.exp(-frequencies*t), 0, np.inf,
                epsabs=cfg.tol, epsrel=cfg.rel_tol, norm="max", limit=1000)[0]),
        }
        if name.startswith("diffusion"):
            methods["SciPy quad_vec t=u^2"] = lambda: Output(integrate.quad_vec(
                lambda u: 2*np.exp(-frequencies*u*u), 0, np.inf,
                epsabs=cfg.tol, epsrel=cfg.rel_tol, norm="max", limit=1000)[0])
        results[name] = {"frequencies": [[s.real, s.imag] for s in frequencies],
                         "methods": measure(methods, truth, repeats)}
        print(f"forward {name} finished", flush=True)
    for count in (2000, 10000):
        t = np.linspace(0, 20, count)
        f = np.exp(-t)
        frequencies = np.array([.1, 1, 3+1j, 10+20j])
        truth = -np.expm1(-(frequencies+1)*t[-1])/(frequencies+1)

        def sampled(method):
            integrand = f[:, None]*np.exp(-t[:, None]*frequencies)
            return Output(method(integrand, x=t, axis=0))

        trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
        results[f"sampled_exp_{count}"] = {"samples": len(t), "methods": measure({
            "LaplaceX spline": lambda: Output(lx.laplace(f, frequencies, t=t)),
            "SciPy simpson": lambda: sampled(integrate.simpson),
            "NumPy trapezoid": lambda: sampled(trapezoid),
        }, truth, repeats)}
    return results


def inverse_cases():
    return [
        ("exp", lambda s: 1/(s+2), lambda t: mp.exp(-2*t), ([1], [1, 2])),
        ("repeated_pole", lambda s: 1/(s+1)**2, lambda t: t*mp.exp(-t),
         ([1], [1, 2, 1])),
        ("damped_cos", lambda s: (s+.5)/((s+.5)**2+9),
         lambda t: mp.exp(-t/2)*mp.cos(3*t), ([1, .5], [1, 1, 9.25])),
        ("sin10", lambda s: 10/(s*s+100), lambda t: mp.sin(10*t),
         ([10], [1, 0, 100])),
        ("diffusion", lambda s: s**(-.5), lambda t: 1/mp.sqrt(mp.pi*t), None),
        ("delayed_step", lambda s: (mp.exp(-s) if isinstance(s, (mp.mpf, mp.mpc))
                                    else cmath.exp(-s))/s,
         lambda t: 0 if t < 1 else 1, None),
    ]


def inverse_suite(repeats, points):
    results = {}
    t = np.linspace(.1, 5, points)
    for name, image, exact, transfer_function in inverse_cases():
        truth = reference(exact, t)

        def mp_inverse(method="cohen", degree=None, dps=15):
            with mp.workdps(dps):
                options = {"method": method}
                if degree is not None:
                    options["degree"] = degree
                values = [float(mp.invertlaplace(image, float(ti), **options)) for ti in t]
            return Output(np.array(values))

        methods = {
            "LaplaceX N64": lambda: Output(lx.inverse(image, t, N=64)),
            "LaplaceX N256": lambda: Output(lx.inverse(image, t, N=256)),
            "mpmath Cohen degree32": lambda: mp_inverse("cohen", 32),
            "mpmath Cohen degree64": lambda: mp_inverse("cohen", 64),
            "mpmath Talbot degree32": lambda: mp_inverse("talbot", 32),
            "mpmath deHoog degree32": lambda: mp_inverse("dehoog", 32),
            "mpmath Stehfest degree32": lambda: mp_inverse("stehfest", 32),
            "mpmath Cohen dps30 auto": lambda: mp_inverse("cohen", None, 30),
        }
        if transfer_function is not None:
            methods["SciPy impulse (rational coefficients)"] = lambda: Output(
                signal.impulse(transfer_function, T=t)[1])
        results[name] = {"times": t.tolist(), "methods": measure(methods, truth, repeats)}
        print(f"inverse {name} finished", flush=True)
    return results


def fine_config():
    cfg = lx.DecomposeConfig()
    cfg.scale_min, cfg.scale_max = .05, .5
    cfg.n_scales, cfg.min_ridge_len = 32, 3
    return cfg


def scipy_fit(t, f, known_count=None):
    """Positive exp sums with analytic Jacobian, deterministic multistart.

    Automatic mode uses the same 1e-6 relative residual-energy stopping rule
    as LaplaceX's fallback, with a maximum of five components. Oracle mode
    receives ONLY the true component count, never true rates or amplitudes.
    """
    amplitude_scale = max(abs(f))
    y = f/amplitude_scale
    lo, hi = .25/t[-1], 4/t[0]
    best = None
    iterations = 0
    for count in ([known_count] if known_count is not None else range(1, 6)):
        model_best = None
        for margin in (1., 4.):
            rates = np.geomspace(margin/t[-1], 1/(margin*t[0]), count)
            rates = np.clip(rates, lo*1.001, hi/1.001)
            start = np.r_[np.full(count, math.log(max(y[0], .01)/count)), np.log(rates)]

            def design(parameters):
                amps, rates = np.exp(parameters[:count]), np.exp(parameters[count:])
                return amps*np.exp(-t[:, None]*rates), rates

            def residual(parameters):
                values, _ = design(parameters)
                return values.sum(axis=1)-y

            def jacobian(parameters):
                values, rates = design(parameters)
                return np.column_stack((values, -t[:, None]*rates*values))

            result = optimize.least_squares(
                residual, start, jac=jacobian, method="trf",
                bounds=(np.r_[np.full(count, -32.), np.full(count, math.log(lo))],
                        np.r_[np.full(count, 8.), np.full(count, math.log(hi))]),
                ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=500,
            )
            iterations += result.nfev
            rss = float(np.dot(result.fun, result.fun))
            if model_best is None or rss < model_best[0]:
                values, fitted_rates = design(result.x)
                model_best = rss, values.sum(axis=1)*amplitude_scale, fitted_rates, result.success
        if best is None or model_best[0] < best[0]:
            best = model_best
        if model_best[0] <= 1e-6*np.dot(y, y):
            break
    return Output(best[1], best[2], {"optimizer_success": bool(best[3]),
                                    "total_nfev": iterations})


def decomposition_cases():
    specs = [
        ("single_uniform", "linear", .01, 15, [2], [1.5], 0, 0),
        ("two_uniform", "linear", .01, 15, [.3, 3], [2, .5], 0, 0),
        ("three_log_pilot", "log", .001, 1000, [.01, 1, 100], [1, 1, 1], 0, 0),
        ("two_log_pilot", "log", .001, 1000, [.01, 100], [1, 1], 0, 0),
        ("three_log_holdout", "log", .002, 400, [.03, 2, 60], [.7, 1.3, .5], 0, 0),
        ("close_rates_holdout", "linear", .01, 20, [.5, .7], [2, .5], 0, 0),
    ]
    for seed in (7, 19, 41):
        specs.append((f"noisy_log_seed{seed}", "log", .002, 400,
                      [.03, 2, 60], [.7, 1.3, .5], .001, seed))
    for name, grid, low, high, rates, amps, sigma, seed in specs:
        t = (np.geomspace if grid == "log" else np.linspace)(low, high, 600)
        clean = sum(a*np.exp(-r*t) for a, r in zip(amps, rates))
        noisy = clean+np.random.default_rng(seed).normal(0, sigma, len(t))
        yield name, t, clean, noisy, np.array(rates), amps, sigma, seed


def decomposition_suite(repeats):
    results = {}
    for name, t, clean, observed, true_rates, amplitudes, sigma, seed in decomposition_cases():
        default, fine = lx.DecomposeConfig(), fine_config()

        def components(output):
            if not isinstance(output, lx.ExponentialComponents):
                output = lx.ExponentialComponents(output)
            return Output(output.reconstruct(t), output.rates)

        methods = {
            "LaplaceX default": lambda: components(lx.decompose(t, observed)),
            "LaplaceX fine wavelets": lambda: components(lx.decompose(t, observed, cfg=fine)),
            "LaplaceX fine wavelets only": lambda: components(lx.decompose(
                t, observed, cfg=fine, auto_deflate=False)),
            "LaplaceX no wavelets": lambda: components(native.decompose_deflate(
                t, observed, default, 5, 1e-6)),
            "SciPy least_squares auto K": lambda: scipy_fit(t, observed),
            "SciPy least_squares known K": lambda: scipy_fit(t, observed, len(true_rates)),
        }
        results[name] = {
            "time_min": float(t[0]), "time_max": float(t[-1]), "samples": len(t),
            "rates": true_rates.tolist(), "amplitudes": list(amplitudes),
            "noise_sigma": sigma, "noise_seed": seed,
            "methods": measure(methods, clean, repeats, observed, true_rates),
        }
        print(f"decomposition {name} finished", flush=True)
    return results


def metadata(args):
    return {
        "python": sys.version, "platform": platform.platform(),
        "processor": platform.processor(), "numpy": np.__version__,
        "scipy": scipy.__version__, "mpmath": mp.__version__,
        "laplacex": lx.__version__, "laplacex_module": lx.__file__,
        "extension_sha256": hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "repeats": args.repeats, "points": args.points,
        "timing": "one warmup, randomized round robin, median; setup included; one BLAS thread",
        "reference": "analytic closed forms, inverse evaluated at 60 dps then rounded to float64",
        "transform_target": "all points abs(error) <= 1e-10 + 1e-8*abs(reference)",
        "fine_wavelet_config": {"scale_min": .05, "scale_max": .5,
                                "n_scales": 32, "min_ridge_len": 3},
        "ablation": "native.decompose_deflate bypasses CWT but uses the same compiled QR/ALS fallback",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--section", choices=["all", "forward", "inverse", "decomposition"],
                        default="all")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--points", type=int, default=64)
    parser.add_argument("--output-dir", type=Path, default=HERE/"comparison")
    args = parser.parse_args()
    if args.repeats < 1 or args.points < 2:
        parser.error("repeats must be positive and points must be at least two")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for section in (["forward", "inverse", "decomposition"]
                    if args.section == "all" else [args.section]):
        results = (forward_suite(args.repeats, args.points) if section == "forward"
                   else inverse_suite(args.repeats, args.points) if section == "inverse"
                   else decomposition_suite(args.repeats))
        report = {"metadata": metadata(args), "results": results}
        (args.output_dir/f"{section}.json").write_text(
            json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        print(f"Saved {section}.json", flush=True)
