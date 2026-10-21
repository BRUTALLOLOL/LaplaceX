"""Generate comparison tables and a static figure from saved measurements."""

from __future__ import annotations

import json
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT/"comparison"


def load(section):
    return json.loads((DATA/f"{section}.json").read_text(encoding="utf-8"))


def successful(methods):
    return {name: row for name, row in methods.items() if row["status"] == "ok"}


def fastest_accurate(methods, prefix):
    eligible = [(name, row) for name, row in successful(methods).items()
                if name.startswith(prefix) and row["passes_transform_target"]]
    return min(eligible, key=lambda item: item[1]["median_seconds"]) if eligible else None


def inverse_summary(inverse):
    points = inverse["metadata"]["points"]
    lines = [f"| Inverse workload ({points} times) | LaplaceX, ms / max error | Fastest passing mpmath, ms / max error | Speed ratio |",
             "|---|---:|---:|---:|"]
    for name, case in inverse["results"].items():
        left = fastest_accurate(case["methods"], "LaplaceX")
        right = fastest_accurate(case["methods"], "mpmath")
        def cell(value):
            if value is None:
                return "No tested setting passed"
            method, row = value
            return f"{row['median_seconds']*1000:.2f} / {row['max_abs_error']:.2e} ({method})"
        speed = (f"{right[1]['median_seconds']/left[1]['median_seconds']:.1f}×"
                 if left and right else "—")
        lines.append(f"| {name} | {cell(left)} | {cell(right)} | {speed} |")
    return lines


def full_table(results, decomposition=False):
    lines = ["| Case | Method | Median ms [min, max] | Max abs error | "
             + ("Clean RMSE | Modes | Rate error |" if decomposition else "Target passed |"),
             "|---|---|---:|---:|" + ("---:|---:|---:|" if decomposition else "---|")]
    for name, case in results.items():
        for method, row in case["methods"].items():
            if row["status"] != "ok":
                suffix = "— | — | — |" if decomposition else "No |"
                message = row["error"].replace("|", "/")
                lines.append(f"| {name} | {method} | Failed: {message} | — | {suffix}")
                continue
            timing = (f"{1000*row['median_seconds']:.3f} "
                      f"[{1000*row['min_seconds']:.3f}, {1000*row['max_seconds']:.3f}]")
            if decomposition:
                rate_error = row["max_rate_relative_error"]
                rate_text = "Count mismatch" if rate_error is None else f"{rate_error:.2e}"
                suffix = f"{row['rmse']:.2e} | {row['component_count']} | {rate_text} |"
            else:
                suffix = "Yes |" if row["passes_transform_target"] else "No |"
            lines.append(f"| {name} | {method} | {timing} | {row['max_abs_error']:.2e} | {suffix}")
    return lines


def figures(inverse, decomposition):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), layout="constrained")
    labels, speed = [], []
    for name, case in inverse["results"].items():
        left = fastest_accurate(case["methods"], "LaplaceX")
        right = fastest_accurate(case["methods"], "mpmath")
        if left and right:
            labels.append(name.replace("_", " "))
            speed.append(right[1]["median_seconds"]/left[1]["median_seconds"])
    ax = axes[0]
    positions = np.arange(len(labels))
    ax.hlines(positions, 1, speed, color="#aac0e1", linewidth=3)
    ax.scatter(speed, positions, color="#245bb2", s=65, zorder=3)
    ax.set_yticks(positions, labels)
    ax.axvline(1, color="#222222", linestyle="--", linewidth=1)
    for position, ratio in zip(positions, speed):
        ax.text(ratio*1.08, position, f"{ratio:.1f}×", va="center")
    ax.set_xscale("log")
    ax.set_xlim(.8, max(speed)*1.6)
    ax.set_title("Batch inverse: speed at a shared accuracy target")
    ax.set_xlabel("Fastest passing mpmath time / fastest passing LaplaceX time")
    ax.text(0, -.23, f"{inverse['metadata']['points']} times per call; all measured methods retained in the report.\n"
            "Unpassed workloads are omitted here, not counted as speed wins.",
            transform=ax.transAxes, fontsize=9, va="top")

    cases = ["single_uniform", "two_uniform", "three_log_pilot", "two_log_pilot",
             "three_log_holdout", "close_rates_holdout"]
    labels = ["single\nuniform", "two\nuniform", "three log\npilot", "two log\npilot",
              "three log\nholdout", "close rates\nholdout", "noisy log\nmean of 3 seeds"]
    ratios, wrong_count = [], []
    for name in cases:
        methods = decomposition["results"][name]["methods"]
        fine, bypass = methods["LaplaceX fine wavelets"], methods["LaplaceX no wavelets"]
        ratios.append(bypass["median_seconds"]/fine["median_seconds"])
        wrong_count.append(not fine["correct_component_count"])
    noisy = [case["methods"] for name, case in decomposition["results"].items()
             if name.startswith("noisy_")]
    ratios.append(float(np.mean([m["LaplaceX no wavelets"]["median_seconds"]
                                /m["LaplaceX fine wavelets"]["median_seconds"] for m in noisy])))
    wrong_count.append(any(not m["LaplaceX fine wavelets"]["correct_component_count"] for m in noisy))
    ax = axes[1]
    bars = ax.bar(np.arange(len(labels)), ratios,
                  color=["#bd5656" if wrong else "#268675" for wrong in wrong_count])
    ax.bar_label(bars, labels=[f"{x:.2f}×" for x in ratios], padding=4, fontsize=9)
    ax.axhline(1, color="#222222", linestyle="--", linewidth=1)
    ax.set_xticks(np.arange(len(labels)), labels, fontsize=8)
    ax.set_ylim(0, max(1.45, max(ratios)*1.2))
    ax.set_ylabel("No-wavelet fitting time / fine-wavelet fitting time")
    ax.set_title("Wavelet ablation: local gains, clear counterexamples")
    ax.text(0, -.23, "Above 1 = wavelets faster; below 1 = wavelets slower.\n"
            "Red = fine-wavelet model returned the wrong component count.",
            transform=ax.transAxes, fontsize=9, va="top")
    fig.savefig(DATA/"comparison.png", dpi=170, bbox_inches="tight")
    fig.savefig(DATA/"comparison.svg", bbox_inches="tight")
    plt.close(fig)


def write_report(forward, inverse, decomposition):
    meta = inverse["metadata"]
    lines = ["# Comparison with Python alternatives", "",
             "Measured locally on 20 October 2026. Every published time is the median of",
             f"{meta['repeats']} measured runs after one warmup; raw run times, failures,",
             "warnings, input definitions, versions, and hashes are preserved in JSON.", "",
             "![Comparison and wavelet ablation](comparison/comparison.png)", "",
             "## Findings", "",
             "LaplaceX has a workload-specific speed advantage over mpmath for batches",
             "of smooth real inverse-transform values. The advantage comes from a compiled",
             "Weeks expansion whose coefficients are reused across times. That inverse",
             "path does not perform wavelet analysis. This is not evidence that wavelets",
             "are responsible for the inverse-transform speed advantage.", "",
             *inverse_summary(inverse), "",
             "The speed ratio compares the fastest **passing** tested configuration in",
             "each package, selected separately for each case. Passing means every point",
             "satisfies `abs(error) <= 1e-10 + 1e-8*abs(reference)`. No speed ratio is",
             "claimed if either package has no passing configuration. Full results below",
             "include all tested configurations and failures, including difficult cases.", "",
             "SciPy is competitive or faster on several forward and fitting workloads.",
             "A general claim of speed or accuracy supremacy is not supported.", "",
             "## What the wavelet ablation establishes", "",
             "The ablation bypasses CWT by calling the existing native `decompose_deflate`",
             "entry point directly. It keeps the same compiled QR/ALS fitting code, sample",
             "arrays, residual stopping tolerance, and AIC setting. Fallback fitting is",
             "limited to five components; CWT determines its own initial mode count.",
             "The no-wavelet route starts its own rate grid search; it receives no true",
             "rates or amplitudes. This tests the value of wavelet initialization relative",
             "to the existing fallback, not an unrelated implementation language. It removes",
             "the log-CWT stage as a whole; it does not separately isolate logarithmic",
             "coordinates from the wavelet operation.", "",
             "Default CWT often finds only one component on mixtures; the public API",
             "then repeats the no-wavelet fallback. Identical results with extra CWT",
             "cost are evidence against attributing those fits to wavelet analysis.", "",
             "The finer preset (`scale_min=0.05`, `scale_max=0.5`, `n_scales=32`,",
             "`min_ridge_len=3`) was explored on the `two_uniform`, `three_log_pilot`,",
             "and `two_log_pilot` signals before the final run. It is an exploratory",
             "configuration, not a universally tuned default. Additional rates, a close",
             "pair, and three noise seeds test its transfer beyond those pilots.", "",
             "On some three-component log-grid cases the fine-wavelet route avoids the",
             "independent rate-search fallback and modestly reduces fitting time. The",
             "counterexamples matter: single- and two-component cases can be slower or",
             "over-segmented, and noisy-fit errors can be effectively identical with or",
             "without wavelets. SciPy fitting remains faster on the tested mixtures.", "",
             "Logarithmic wavelets are therefore a configurable initialization mechanism",
             "with local measured gains, not a demonstrated general accuracy advantage.", "",
             "## Fairness and precision", "",
             f"- Python: `{meta['python'].splitlines()[0]}`; platform: `{meta['platform']}`.",
             f"- NumPy `{meta['numpy']}`, SciPy `{meta['scipy']}`, mpmath `{meta['mpmath']}`.",
             "- Benchmark sections run serially; each uses one BLAS thread. Method order",
             "  is shuffled deterministically on each repetition. Setup and array",
             "  conversion inside each method are timed; imports and truth generation are not.",
             "- LaplaceX uses the compiled Release wheel from the preceding fixes:",
             "  Zig 0.16.0 / Clang 21.1.0, C++17, `-O3`, Windows GNU x64 target.",
             "  The extension SHA-256 and benchmark-script hash are recorded in each JSON.",
             "- Inverse references use analytic formulas evaluated at 60 decimal digits,",
             "  then rounded to float64. All outputs are compared as float64. This",
             "  measures float64 workloads; mpmath can deliver much higher precision.",
             "  A displayed zero error means agreement after float64 rounding, not an",
             "  exact numerical inversion.",
             "- mpmath includes Cohen at degrees 32 and 64, Talbot/de Hoog/Stehfest at",
             "  degree 32, and automatic Cohen at requested 30 dps. Fixed-degree calls",
             "  begin at 15 dps but mpmath raises its working precision internally.",
             "  Equal term counts do not imply equal work or equal arithmetic precision.",
             f"  The local mpmath backend is `{mp.libmp.BACKEND}`; optional gmpy2 acceleration",
             "  was not installed. Results do not compare every possible package build.",
             "- A separate follow-up adds LaplaceX N512 on every inverse case and de Hoog",
             "  degree64 on the delayed step after the initial settings exposed failures.",
             "  Original runs are retained; follow-up settings are labeled and saved separately.",
             "- Both LaplaceX and mpmath inverse calls receive only a callable image.",
             "  `scipy.signal.impulse` receives numerator/denominator coefficients and",
             "  is explicitly labeled as a richer-input rational-system comparison.",
             "- Forward callable comparisons use shared absolute/relative tolerances",
             "  `1e-10`/`1e-8`. SciPy integrates to infinity; LaplaceX truncates at `1e6`.",
             "  Omitted tails are negligible for these specific frequencies and signals.",
             "  The diffusion test also gives SciPy the same squared-time endpoint",
             "  regularization used by LaplaceX, as an explicitly labeled variant.",
             "- Sampled-data methods receive the same samples and finite-window truth.",
             "  They are evaluated for speed and approximation error, not treated as",
             "  interchangeable algorithms with identical accuracy guarantees.",
             "- Decomposition judges predictions against the clean signal even for",
             "  noisy inputs; it also reports true/recovered counts and rate errors.",
             "  SciPy fitting uses positive amplitudes, an analytic Jacobian, bounded",
             "  log rates, two deterministic starts, and tight optimizer tolerances.",
             "  Automatic K uses the same residual-energy stopping rule and five-mode",
             "  budget. Known-K receives the true count only, not rates or amplitudes.",
             "- These synthetic cases and three noise seeds are a bounded local study,",
             "  not confidence intervals for all inputs, CPUs, or real experimental data.", "",
             "## All forward measurements", "", *full_table(forward["results"]), "",
             "## All inverse measurements", "", *full_table(inverse["results"]), "",
             "## All decomposition measurements", "",
             "Rate error is the largest relative difference after sorting true and fitted",
             "rates; it is omitted when component counts disagree. Clean RMSE is computed",
             "against the generating signal. `pilot` labels identify explored cases.", "",
             *full_table(decomposition["results"], True), "",
             "## Reproduce", "", "```bash", 'python -m pip install ".[benchmark]"',
             "python benchmarks/compare.py --section all --repeats 5 --points 64",
             "python benchmarks/compare_extended.py --repeats 5 --points 64",
             "python benchmarks/render_comparison.py", "```", "",
             "Raw data: [forward](comparison/forward.json), [inverse](comparison/inverse.json),",
             "[inverse follow-up](comparison/inverse_extended.json),",
             "[decomposition](comparison/decomposition.json). The script does not change",
             "library code or automatically select favorable cases. Existing output files",
             "are replaced when the corresponding section is rerun.", "",
             "## Primary documentation", "",
             "- [mpmath inverse algorithms and internal precision](https://mpmath.org/doc/current/calculus/inverselaplace.html)",
             "- [SciPy quad tolerances and complex integrands](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html)",
             "- [SciPy quad_vec](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad_vec.html)",
             "- [SciPy least_squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)",
             "- [SciPy rational-system impulse response](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.impulse.html)", ""]
    (ROOT/"COMPARISON.md").write_text("\n".join(lines), encoding="utf-8")


def update_readme(inverse, decomposition):
    lines = ["<!-- comparison:start -->", "## Comparison with Python alternatives", "",
             "Local measurements compare LaplaceX with mpmath, SciPy, and NumPy.",
             "The strongest result is fast **batch inversion of smooth signals**:",
             "one compiled Weeks expansion serves many evaluation times. Every speed",
             "ratio below compares the fastest tested configuration in each package",
             "that met the same pointwise accuracy target,",
             "`abs(error) <= 1e-10 + 1e-8*abs(reference)`.", "",
             *inverse_summary(inverse), "",
             "Times are medians of five runs after warmup on Windows x64 / Python 3.14;",
             "errors are maximum absolute differences from 60-digit analytic references",
             "rounded to float64. mpmath uses higher internal precision and can provide",
             "accuracy beyond float64. The local mpmath uses its Python backend; optional",
             "gmpy2 acceleration was not installed. All tested settings and failures are published.", "",
             "**Wavelet-specific evidence.** Removing the logarithmic CWT stage and",
             "calling the same compiled iterative fallback gives this controlled ablation:", "",
             "| Three-component log-grid signal | Fine log-wavelets, ms | No wavelets, ms | Speed ratio |",
             "|---|---:|---:|---:|"]
    for case_name, label in [("three_log_pilot", "Rates 0.01, 1, 100 (pilot)"),
                             ("three_log_holdout", "Rates 0.03, 2, 60 (holdout)")]:
        methods = decomposition["results"][case_name]["methods"]
        fine, bypass = methods["LaplaceX fine wavelets"], methods["LaplaceX no wavelets"]
        ratio = bypass["median_seconds"]/fine["median_seconds"]
        lines.append(f"| {label} | {fine['median_seconds']*1000:.2f} | "
                     f"{bypass['median_seconds']*1000:.2f} | {ratio:.2f}× |")
    lines += ["",
              "The fine preset uses `scale_min=0.05`, `scale_max=0.5`, `n_scales=32`,",
              "and `min_ridge_len=3`; it was explored on pilot signals before testing",
              "additional rates and noise seeds. It can avoid the independent grid-search",
              "fallback on some mixtures. This is a **local speed benefit**, not a general",
              "accuracy advantage: noisy-fit errors were effectively unchanged, and the",
              "same preset overestimated the mode count on other signals.", "",
              "![Batch inverse comparison and wavelet ablation](benchmarks/comparison/comparison.png)", "",
              "The counterexamples are part of the result: SciPy was faster on the larger",
              "forward batch and on the tested exponential fits. Its `impulse` solver was",
              "also faster when given rational transfer-function coefficients. Default",
              "LaplaceX wavelets often defer to the non-wavelet fallback, and no tested",
              "LaplaceX inverse setting met the target for the diffusion or delayed-step",
              "examples. The Weeks inverse and forward quadrature paths do not use wavelets,",
              "so their speed gains cannot be attributed to logarithmic wavelet analysis.", "",
              "[Full comparison, settings, limitations, and raw results](benchmarks/COMPARISON.md)",
              "are reproducible with:", "", "```bash", 'python -m pip install ".[benchmark]"',
              "python benchmarks/compare.py --section all --repeats 5 --points 64",
              "python benchmarks/compare_extended.py --repeats 5 --points 64",
              "python benchmarks/render_comparison.py --update-readme", "```", "",
              "<!-- comparison:end -->", ""]
    path = ROOT.parent/"README.md"
    text = path.read_text(encoding="utf-8")
    start_marker, end_marker = "<!-- comparison:start -->", "<!-- comparison:end -->"
    if start_marker in text:
        start = text.index(start_marker)
        end = text.index(end_marker, start)+len(end_marker)
        text = text[:start]+"\n".join(lines).rstrip()+text[end:]
    else:
        target = "## Quick Start"
        text = text.replace(target, "\n".join(lines)+"\n"+target, 1)
    path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--update-readme", action="store_true")
    args = parser.parse_args()
    forward, inverse, decomposition = [load(section) for section in
                                       ("forward", "inverse", "decomposition")]
    extended = load("inverse_extended")
    for dataset in (forward, decomposition, extended):
        for key in ("extension_sha256", "points", "repeats", "script_sha256"):
            if dataset["metadata"][key] != inverse["metadata"][key]:
                raise ValueError(f"Mismatched benchmark metadata: {key}")
    for name, case in extended["results"].items():
        inverse["results"][name]["methods"].update(case["methods"])
    figures(inverse, decomposition)
    write_report(forward, inverse, decomposition)
    if args.update_readme:
        update_readme(inverse, decomposition)
    print("Saved COMPARISON.md and comparison.png/.svg")
