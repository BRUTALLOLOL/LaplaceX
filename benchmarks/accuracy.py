"""Deterministic analytic accuracy and timing checks for an installed LaplaceX.

Run with the same interpreter, hardware, and compiler for before/after comparisons:
    python benchmarks/accuracy.py --output benchmarks/after.json
Optional independent 60-digit references: pip install mpmath
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

import laplacex as lx

try:
    import mpmath as mp
except ImportError:
    mp = None


def run():
    results = {}

    def measure(name, operation, expected):
        start = perf_counter()
        actual = np.asarray(operation())
        elapsed = perf_counter() - start
        difference = abs(actual - expected)
        results[name] = {
            "max_abs_error": float(np.max(difference)),
            "normalized_max_error": float(np.max(difference) / np.max(abs(expected))),
            "seconds": elapsed,
        }

    frequencies = np.array([.1, 1, 3+1j, 10+20j])
    calls = 0

    def decay(t):
        nonlocal calls
        calls += 1
        return math.exp(-2*t)

    measure("forward_exp", lambda: lx.laplace(decay, frequencies), 1/(frequencies+2))
    results["forward_exp"]["callback_evaluations"] = calls
    cfg = lx.TransformConfig()
    cfg.tol = 1e-20
    measure("forward_fast_decay", lambda: lx.laplace(lambda t: math.exp(-1e12*t),
                                                  1e12, cfg=cfg), 5e-13)
    measure("forward_diffusion", lambda: lx.laplace(lambda t: 1/math.sqrt(t), 2),
            math.sqrt(math.pi/2))
    t = np.linspace(.001, 20, 2000)
    measure("sampled_exp", lambda: lx.laplace(np.exp(-t), frequencies, t=t),
            1/(frequencies+1))
    t = np.linspace(.1, 5, 100)
    # Evaluate closed forms at high precision when mpmath is available.
    if mp:
        mp.mp.dps = 60
        reference = lambda fn: np.array([float(fn(mp.mpf(float(x)))) for x in t])
    else:
        reference = lambda fn: np.array([float(fn(x)) for x in t])
    exp = mp.exp if mp else math.exp
    sin = mp.sin if mp else math.sin
    cos = mp.cos if mp else math.cos
    for name, image, exact in [
        ("inverse_exp", lambda s: 1/(s+2), lambda t: exp(-2*t)),
        ("inverse_step", lambda s: 1/s, lambda t: 1),
        ("inverse_ramp", lambda s: 1/s**2, lambda t: t),
        ("inverse_damped_cos", lambda s: (s+.5)/((s+.5)**2+9),
         lambda t: exp(-.5*t)*cos(3*t)),
    ]:
        measure(name, lambda image=image: lx.inverse(image, t), reference(exact))
    measure("inverse_sine_N256", lambda: lx.inverse(lambda s: 10/(s*s+100), t, N=256),
            reference(lambda t: sin(10*t)))
    measure("inverse_sine_N64", lambda: lx.inverse(lambda s: 10/(s*s+100), t),
            reference(lambda t: sin(10*t)))
    c = lx.ExponentialComponent()
    c.amplitude, c.rate, c.frequency, c.phase = 2, .5, 3, .7
    comps = lx.ExponentialComponents([c])
    s = frequencies
    expected = 2*((s+.5)*math.cos(.7)-3*math.sin(.7))/((s+.5)**2+9)
    measure("oscillatory_component_image", lambda: comps.reconstruct_laplace(s), expected)
    t = np.linspace(.01, 15, 600)
    for name, signal in [
        ("decompose_single", 1.5*np.exp(-2*t)),
        ("decompose_two", 2*np.exp(-.3*t)+.5*np.exp(-3*t)),
    ]:
        measure(name, lambda signal=signal: lx.decompose(t, signal).reconstruct(t), signal)
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "module": lx.__file__,
        "reference_precision_digits": 60 if mp else 16,
        "results": results,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = json.dumps(run(), indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report + "\n", encoding="utf-8")
    print(report)
