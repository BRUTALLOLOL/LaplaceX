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
        fine, bypass = methods["LaplaceX default"], methods["LaplaceX no wavelets"]
        ratios.append(bypass["median_seconds"]/fine["median_seconds"])
        wrong_count.append(not fine["correct_component_count"])
    noisy = [case["methods"] for name, case in decomposition["results"].items()
             if name.startswith("noisy_")]
    ratios.append(float(np.mean([m["LaplaceX no wavelets"]["median_seconds"]
                                /m["LaplaceX default"]["median_seconds"] for m in noisy])))
    wrong_count.append(any(not m["LaplaceX default"]["correct_component_count"] for m in noisy))
    ax = axes[1]
    bars = ax.bar(np.arange(len(labels)), ratios,
                  color=["#bd5656" if wrong else "#268675" for wrong in wrong_count])
    ax.bar_label(bars, labels=[f"{x:.2f}×" for x in ratios], padding=4, fontsize=9)
    ax.axhline(1, color="#222222", linestyle="--", linewidth=1)
    ax.set_xticks(np.arange(len(labels)), labels, fontsize=8)
    ax.set_ylim(0, max(1.45, max(ratios)*1.2))
    ax.set_ylabel("No-wavelet fitting time / default log-wavelet fitting time")
    ax.set_title("Wavelet ablation: same optimizer, different initialization")
    ax.text(0, -.23, "Above 1 = wavelets faster; below 1 = wavelets slower.\n"
            "Red = default wavelet model returned the wrong component count.",
            transform=ax.transAxes, fontsize=9, va="top")
    fig.savefig(DATA/"comparison.png", dpi=170, bbox_inches="tight")
    fig.savefig(DATA/"comparison.svg", bbox_inches="tight")
    plt.close(fig)


def decomposition_pass(case, row):
    if row["status"] != "ok" or not row["correct_component_count"]:
        return False
    if case["noise_sigma"]:
        return row["rmse"] <= case["noise_sigma"] and row["max_rate_relative_error"] <= .1
    return row["normalized_max_error"] <= 1e-8 and row["max_rate_relative_error"] <= 1e-5


def best_scipy(case, known=False):
    candidates = [(name, row) for name, row in case["methods"].items()
                  if name.startswith("SciPy") and
                  (("known K" in name) if known else ("auto K" in name))
                  and decomposition_pass(case, row)]
    return min(candidates, key=lambda pair: pair[1]["median_seconds"]) if candidates else None


def decomposition_summary(decomposition):
    lines = ["| Signal | LaplaceX default, ms | Fastest passing SciPy auto-K, ms | Speed ratio | Modes: LaplaceX / truth |",
             "|---|---:|---:|---:|---:|"]
    for name, case in decomposition["results"].items():
        ours = case["methods"]["LaplaceX default"]
        scipy = best_scipy(case)
        other = f"{1000*scipy[1]['median_seconds']:.2f}" if scipy else "No passing configuration"
        ratio = (f"{scipy[1]['median_seconds']/ours['median_seconds']:.2f}×"
                 if scipy and decomposition_pass(case, ours) else "—")
        lines.append(f"| {name} | {1000*ours['median_seconds']:.2f} | {other} | {ratio} | "
                     f"{ours['component_count']} / {len(case['rates'])} |")
    return lines


def weak_summary(decomposition):
    lines = ["| Weak mode | LaplaceX modes / max abs error | SciPy residual-stop modes / max abs error | SciPy BIC / known-K accuracy |",
             "|---|---:|---:|---|"]
    for name in ("weak_fast", "weak_middle", "weak_slow"):
        case = decomposition["results"][name]
        ours = case["methods"]["LaplaceX default"]
        legacy = case["methods"]["SciPy least_squares auto K"]
        strong, oracle = best_scipy(case), best_scipy(case, True)
        strongest_error = min(pair[1]["max_abs_error"] for pair in (strong, oracle) if pair)
        lines.append(f"| {name} | {ours['component_count']} / {ours['max_abs_error']:.2e} | "
                     f"{legacy['component_count']} / {legacy['max_abs_error']:.2e} | "
                     f"Also recovers 3 modes; error down to {strongest_error:.2e} |")
    return lines


def ablation_text(decomposition):
    def ratio(name):
        methods = decomposition["results"][name]["methods"]
        return methods["LaplaceX no wavelets"]["median_seconds"]/methods["LaplaceX default"]["median_seconds"]
    return [
        f"With the same optimizer, default wavelet initialization is **{ratio('three_log_pilot'):.2f}×** faster",
        f"on the widely separated three-mode case and **{ratio('three_log_holdout'):.2f}×** faster on its holdout.",
        f"The two-uniform and close-rate ratios are **{ratio('two_uniform'):.2f}×** and",
        f"**{ratio('close_rates_holdout'):.2f}×** respectively; above 1 means wavelets are faster.",
        "These ratios compare the native rate-grid route",
        "with the public wavelet route, including its Python routing overhead.",
    ]


def pair_optimization_summary(decomposition):
    before = json.loads((ROOT/"comparison-before-pair-optimization/decomposition.json").read_text())
    lines = ["| Previously slower case | Before: wavelet / grid ms | Now: wavelet / grid ms | Wavelet speedup vs previous | Now: grid / wavelet |",
             "|---|---:|---:|---:|---:|"]
    for name in ("two_uniform", "close_rates_holdout"):
        old = before["results"][name]["methods"]
        new = decomposition["results"][name]["methods"]
        ow, og = [old[k]["median_seconds"] for k in ("LaplaceX default", "LaplaceX no wavelets")]
        nw, ng = [new[k]["median_seconds"] for k in ("LaplaceX default", "LaplaceX no wavelets")]
        lines.append(f"| {name} | {1000*ow:.2f} / {1000*og:.2f} | {1000*nw:.2f} / {1000*ng:.2f} | {ow/nw:.2f}× | {ng/nw:.2f}× |")
    return lines


def fresh_validation_text():
    fresh = load("validation_fresh")
    old = json.loads((ROOT/"comparison-before-noisy-fix/validation_fresh.json").read_text())
    holdout = load("validation_holdout")
    primary = load("validation")
    changed = sum(len(a["recovered_rates"]) != len(b["recovered_rates"])
                  for a, b in zip(old["results"], fresh["results"]))
    def noisy_p95(dataset):
        return 100*np.quantile([r["max_rate_relative_error"] for r in dataset["results"]
                                if r["noise_sigma"] and r["correct_count"]], .95)
    worst = max(r["max_rate_relative_error"] for d in (primary, fresh, holdout)
                for r in d["results"] if r["noise_sigma"] and r["correct_count"])
    return [
        f"The two noisy failures in seed {fresh['seed']} are fixed: **{old['correct_counts']}/{old['cases']} → {fresh['correct_counts']}/{fresh['cases']}** correct counts.",
        f"There are **{changed} changed count decisions** on that regression set.",
        f"An additional seed ({holdout['seed']}) recovers **{holdout['correct_counts']}/{holdout['cases']}** counts.",
        f"Across the three seeds: **{sum(d['correct_counts'] for d in (primary, fresh, holdout))}/{sum(d['cases'] for d in (primary, fresh, holdout))}**, with fallback disabled.",
        f"The formerly failing set's noisy rate-error 95th percentile is **{noisy_p95(old):.2f}% before / {noisy_p95(fresh):.2f}% after**.",
        f"The additional seed's noisy rate-error 95th percentile is **{noisy_p95(holdout):.2f}%**.",
        f"Across all noisy mixtures, the largest rate error is **{100*worst:.2f}%**;",
        "this pre-existing limitation shows why correct counts do not guarantee precise individual rates.",
        "Seed 20261107 is now regression data because its failures informed these fixes;",
        "seed 20261119 was chosen before implementation and was not used to tune the rules.",
    ]


def noisy_fix_summary():
    before = json.loads((ROOT/"comparison-before-noisy-fix/validation_fresh.json").read_text())
    after = load("validation_fresh")
    lines = ["| Previously failing case | Modes before → after / truth | Clean RMSE before → after | Largest rate error after |",
             "|---|---:|---:|---:|"]
    for a, b in zip(before["results"], after["results"]):
        if a["correct_count"]:
            continue
        assert b["correct_count"]
        lines.append(f"| {a['grid']} grid, case {a['index']} | {len(a['recovered_rates'])} → {len(b['recovered_rates'])} / {len(b['rates'])} | "
                     f"{a['clean_rmse']:.2e} → {b['clean_rmse']:.2e} | {100*b['max_rate_relative_error']:.2f}% |")
    return lines


def pair_check_summary():
    pairs = load("pairs")
    lines = [f"A separate **{pairs['metadata']['repeats']}-run** randomized check of the final binary confirms the local pair timings.",
             "Cells show median [25th, 75th percentile] milliseconds; every fit must pass the component-count and accuracy checks.", "",
             "| Case | Default wavelets, ms | Rate grid, ms | Grid / wavelet |",
             "|---|---:|---:|---:|"]
    for name, case in pairs["results"].items():
        for row in case["methods"].values():
            if not decomposition_pass(case, row):
                raise ValueError(f"Pair check failed count/accuracy: {name}")
        wave, grid = [case["methods"][key] for key in ("LaplaceX default", "LaplaceX no wavelets")]
        def cell(row):
            lower, upper = np.quantile(row["timings_seconds"], [.25, .75])*1000
            return f"{1000*row['median_seconds']:.2f} [{lower:.2f}, {upper:.2f}]"
        lines.append(f"| {name} | {cell(wave)} | {cell(grid)} | {grid['median_seconds']/wave['median_seconds']:.2f}× |")
    return lines


def validation_text(validation):
    clean = [r for r in validation["results"] if r["noise_sigma"] == 0]
    noisy = [r for r in validation["results"] if r["noise_sigma"] > 0]
    errors = [r["max_rate_relative_error"] for r in noisy if r["correct_count"]]
    return [
        f"The separately seeded validation recovered the correct count in **{validation['correct_counts']}/{validation['cases']}** mixtures",
        "with the rate-grid fallback disabled. It includes 80 clean and 80 noisy signals,",
        "1–4 modes, linear and logarithmic sampling, and independently varied amplitudes/rates.",
        f"Maximum clean rate error: `{max(r['max_rate_relative_error'] for r in clean if r['correct_count']):.3e}` relative.",
        f"For noise σ=0.001, median / 95th-percentile / maximum rate errors were",
        f"**{100*np.median(errors):.2f}% / {100*np.quantile(errors,.95):.2f}% / {100*max(errors):.2f}%**.",
        "Correct counts do not imply exact rates in noisy, poorly conditioned mixtures.",
        "These are finite synthetic validation results, not a guarantee for all signals.",
    ]


def wavelet_figure(decomposition):
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), layout="constrained")
    labels, ratios = [], []
    for name, case in decomposition["results"].items():
        other = best_scipy(case)
        ours = case["methods"]["LaplaceX default"]
        if other and decomposition_pass(case, ours):
            labels.append(name.replace("_", " "))
            ratios.append(other[1]["median_seconds"]/ours["median_seconds"])
    ax = axes[0]
    y = np.arange(len(labels))
    ax.hlines(y, 1, ratios, color="#bed9d1", linewidth=3)
    ax.scatter(ratios, y, color="#268675", s=40, zorder=3)
    for i, ratio in enumerate(ratios):
        ax.text(ratio*1.08, i, f"{ratio:.1f}×", va="center", fontsize=9)
    ax.set_yticks(y, labels, fontsize=8)
    ax.axvline(1, color="black", linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlim(min(.8, min(ratios)*.8), max(ratios)*1.8)
    ax.set_title("Component recovery: speed at matched accuracy")
    ax.set_xlabel("Fastest passing SciPy automatic-K time / LaplaceX time")
    ax.text(0, -.16, "TRF and LM, including BIC selection; same observations, no true K.\n"
            "Known-K oracle results are reported separately.", transform=ax.transAxes,
            fontsize=9, va="top")
    names = ["weak_fast", "weak_middle", "weak_slow"]
    palette = ["#245bb2", "#bd5656", "#7d8b8e", "#9b83bd"]
    series = [[], [], [], []]
    counts = [[], [], [], []]
    for name in names:
        case = decomposition["results"][name]
        rows = [case["methods"]["LaplaceX default"],
                case["methods"]["SciPy least_squares auto K"],
                best_scipy(case)[1], best_scipy(case, True)[1]]
        for i, row in enumerate(rows):
            series[i].append(max(1e-17, row["max_abs_error"]))
            counts[i].append(row["component_count"])
    ax = axes[1]
    for j, label in enumerate(["LaplaceX default", "SciPy residual-stop auto K",
                               "SciPy best passing auto K", "SciPy best known-K oracle"]):
        positions = np.arange(3)+(j-1.5)*.2
        bars = ax.bar(positions, series[j], width=.18, color=palette[j], label=label)
        ax.bar_label(bars, labels=[f"K={k}" for k in counts[j]], padding=3, fontsize=8)
    ax.set_yscale("log")
    ax.set_ylim(1e-17, .02)
    ax.set_xticks(range(3), ["weak fast", "weak middle", "weak slow"])
    ax.set_ylabel("Maximum absolute reconstruction error")
    ax.set_title("1000:1 amplitudes: missing a mode versus recovering it", pad=60)
    ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.02), ncol=2, fontsize=8)
    ax.text(0, -.16, "Stronger SciPy selection also reaches float64 roundoff.\n"
            "This accuracy gain concerns automatic mode selection, not extra precision.",
            transform=ax.transAxes, fontsize=9, va="top")
    fig.savefig(DATA/"wavelet_comparison.png", dpi=170, bbox_inches="tight")
    fig.savefig(DATA/"wavelet_comparison.svg", bbox_inches="tight")
    plt.close(fig)


def write_report(forward, inverse, decomposition, validation):
    meta = decomposition["metadata"]
    before = json.loads((ROOT/"comparison-before-wavelet/decomposition.json").read_text())
    lines = ["# Log-plane wavelet recovery and Python comparisons", "",
             "Measured on 21 October 2026 after the noisy-component fixes. Times are medians",
             "of five runs after one warmup, with shuffled method order and one BLAS thread.", "",
             "![Wavelet recovery comparison](comparison/wavelet_comparison.png)", "",
             "## Component recovery", "", *decomposition_summary(decomposition), "",
             "A passing automatic-count fit must recover the true count. Clean cases also",
             "require normalized max reconstruction error ≤1e-8 and rate error ≤1e-5.",
             "Noisy cases require clean-signal RMSE ≤noise σ and rate error ≤10%.",
             "The comparator is the fastest passing automatic-K SciPy configuration among",
             "the original residual-stop TRF fitter and TRF/LM BIC fitters. All outputs,",
             "including failures and incorrect counts, remain in the full tables below.", "",
             "## What changed", "",
             "- Analyze `f(exp(x))`, rather than `exp(x)*f(exp(x))`. The negative log",
             "  derivative of one exponential has the profile `A*exp(z-exp(z))`,",
             "  `z=x+log(rate)`: rate changes translate the profile without changing",
             "  its height. The old time factor favored slow components.",
             "- Track persistent Mexican-Hat zero crossings from fine to coarse scales,",
             "  retaining branches when coarse scales merge them. Initialize all persistent",
             "  components in one fit instead of repeatedly fitting incomplete models.",
             "- Fit rates and signed amplitudes jointly with an analytic Jacobian, damped",
             "  Levenberg–Marquardt steps and scaled Householder QR. Final fitting uses",
             "  the original observations, so log-grid interpolation is initialization only.",
             "- Propagate sample noise through interpolation and wavelet filter weights,",
             "  including covariance from reused observations. Condition residual proposals",
             "  on both the existing amplitudes and rates.",
             "- Use BIC by default, refit candidate deletions, and screen nearly dependent",
             "  modes using conditional amplitude uncertainty. AICc remains configurable.",
             "- Respect the component budget in both stages. Handle constants, noise-only",
             "  signals, and large/small amplitude scales without spurious decay modes or",
             "  overflowing the Python residual-energy calculation.", "",
             "- Share an exact single-exponential shortcut between both routes. It accepts",
             "  a model only after fitting every sample at the common numerical floor.",
             "  Its single-mode speed benefit is not attributed to wavelet analysis.", "",
             "## Fixing the remaining noisy count errors", "",
             "A fast component outside the observed characteristic-time window had an",
             "unconstrained extrapolated amplitude despite passing BIC. For uncertain modes",
             "with `rate*t_min > 1`, the existing conditional amplitude test now requires",
             "a refitted deletion. In-window overlapping modes retain the BIC/AICc decision.", "",
             "The missed fourth component was suppressed by the initial interpolation-residual",
             "noise estimate. Once fit residuals are within twice the expected noise energy,",
             "automatic-noise searches cap their proposal noise scale by the residual MAD.",
             "Explicit user-supplied noise is respected; the original estimate is never increased.",
             "The covariance propagation, QR solver and pair-case optimizations remain in place.", "",
             *noisy_fix_summary(), "",
             "## Optimizing the previously slower pair cases", "",
             "- Evaluate and cache exact propagated noise variance only at candidate crossings.",
             "- Pad convolution boundaries once per residual rather than clamp each filter tap.",
             "- For unresolved pairs, search both sides of the strongest residual ridge using",
             "  conditional Jacobian projections. Keep the original proposals and existing selection rules.",
             "- Reduce the Jacobian with Householder QR once per nonlinear iteration, then",
             "  try damping changes on the small triangular problem. No normal equations or",
             "  relaxed tolerances are introduced. The QR improvement is shared by both routes.", "",
             *pair_optimization_summary(decomposition), "",
             *pair_check_summary(), "",
             "Before/after timings come from separate local runs; same-run grid/wavelet ratios",
             "are the cleaner ablation. All models and final fits still use the original observations.", "",
             "## Accuracy of automatic model selection", "", *weak_summary(decomposition), "",
             "The 1000:1 amplitude cases expose the old SciPy automatic baseline's",
             "`1e-6` residual-energy stopping rule: it can accept two modes before a weak",
             "third mode is recovered. LaplaceX finds all three. This is an advantage over",
             "that automatic stopping policy, not universal numerical accuracy superiority:",
             "the stronger SciPy BIC and known-K variants also recover the modes accurately.",
             "Float64 roundoff differences between successful fits are not meaningful accuracy wins.", "",
             "## Validation across three seeds", "", *validation_text(validation), "",
             "[Full validation protocol and individual outcomes](comparison/validation.json).",
             "Model-selection development used seed 20261023; the separate validation",
             "uses seed 20261029. Later sparse-variance and exact-single-mode performance",
             "optimizations were revalidated on the same protocol; no selection parameters",
             "were tuned to the second seed.",
             "Validation timings are single observations and are not used for speed claims.", "",
             *fresh_validation_text(), "",
             "The before-fix results and",
             "[previous critical.cpp](comparison-before-noisy-fix/critical.cpp) are archived",
             "alongside the previously benchmarked extension hash; these are actual previous fits.", "",
             "## Wavelet ablation", "",
             "`LaplaceX no wavelets` bypasses CWT and searches a 64-point log-rate grid.",
             "It shares the new QR/LM solver, noise estimation, conditional candidate test,",
             "BIC and deletion refits. This isolates initialization and candidate generation",
             "from the optimizer rewrite. It removes the entire log-CWT stage; it does not",
             "separately compare linear-time and log-time wavelets.", "",
             "The default public route also includes Python fallback/noise checks. The",
             "native no-wavelet entry omits that routing overhead; the explicitly",
             "wavelet-only rows omit fallback checks and provide an additional comparison.", "",
             *ablation_text(decomposition), "",
             "![Inverse speed and wavelet ablation](comparison/comparison.png)", "",
             "The default configuration now uses 16 scales over [0.05,0.5] and 256 log points.",
             "The older fine preset remains in the matrix. Both automatic fallback behavior",
             "and the explicitly wavelet-only route are tested. Excess residual consistent",
             "with estimated noise no longer triggers a redundant rate-grid fit.", "",
             "## Before and after on identical fixtures", "",
             "These compare the previous default library with the rewritten default library;",
             "they include the combined effects of wavelet recovery, fitting and selection.",
             "The ablation above is the narrower test of wavelet contribution.", "",
             "| Case | Previous default ms / modes | New default ms / modes | Max abs error before → after |",
             "|---|---:|---:|---:|"]
    for name, old_case in before["results"].items():
        new_case = decomposition["results"][name]
        assert old_case["rates"] == new_case["rates"] and old_case["amplitudes"] == new_case["amplitudes"]
        old, new = old_case["methods"]["LaplaceX default"], new_case["methods"]["LaplaceX default"]
        lines.append(f"| {name} | {1000*old['median_seconds']:.2f} / {old['component_count']} | "
                     f"{1000*new['median_seconds']:.2f} / {new['component_count']} | "
                     f"{old['max_abs_error']:.2e} → {new['max_abs_error']:.2e} |")
    lines += ["", "Historical raw files and benchmark scripts are preserved in",
              "[comparison-before-wavelet](comparison-before-wavelet/decomposition.json).", "",
              "## Forward and inverse comparisons", "", *inverse_summary(inverse), "",
              "Inverse ratios select the fastest passing tested setting in each package.",
              "Every point must satisfy `abs(error) <= 1e-10+1e-8*abs(reference)`.",
              "The Weeks inverse and forward quadrature do not use wavelet decomposition.",
              "No wavelet-related accuracy gain is attributed to those paths. The inverse",
              "still has limitations on singular/discontinuous signals. SciPy's rational",
              "impulse response has richer inputs and is retained as a separate comparison.", "",
              "## Fairness, settings, and limits", "",
              f"- Python `{meta['python'].splitlines()[0]}`; `{meta['platform']}`.",
              f"- NumPy `{meta['numpy']}`, SciPy `{meta['scipy']}`, mpmath `{meta['mpmath']}`",
              f"  with its `{mp.libmp.BACKEND}` backend (gmpy2 not installed).",
              "- LaplaceX: Release Windows x64 wheel, Zig 0.16.0 / Clang 21.1.0, C++17, -O3.",
              "  Extension and script SHA-256 hashes, all run times and warnings are in JSON.",
              "- Sections run serially. Setup, optimization and output conversion are timed;",
              "  imports and reference generation are not. Medians are local observations,",
              "  not cross-machine performance guarantees or confidence intervals.",
              "- SciPy fitters use positive amplitudes, analytic Jacobians, two deterministic",
              "  starts, tight tolerances, and at most five components. TRF uses bounded",
              "  log rates/amplitudes; LM uses their unconstrained logarithmic parameterization.",
              "  Known-K variants receive only the true count, never true rates/amplitudes.",
              "- BIC variants evaluate K=1..5, stopping early only when reaching the common",
              "  `1e-28*signal_energy` floor makes larger models incapable of improving BIC.",
              "  The original residual-stop fitter is retained for comparability, not treated",
              "  as the only SciPy alternative. Known-K oracle outcomes are reported below.",
              "- Inverse references are analytic formulas evaluated at 60 digits then rounded",
              "  to float64. Zero displayed error means agreement after rounding. mpmath's",
              "  internal precision can exceed float64, so speed is not equal-precision proof.",
              "- Forward callable tolerances match at 1e-10 absolute and 1e-8 relative.",
              "  LaplaceX truncates at 1e6; omitted tails are negligible here. SciPy also gets",
              "  an explicitly labeled squared-time endpoint variant for diffusion.",
              "- These benchmarks use decaying real exponential mixtures, not arbitrary",
              "  oscillatory signals or relaxation distributions. Decomposition no longer",
              "  fabricates oscillation frequencies from Mexican-Hat side lobes. Manual",
              "  oscillatory component reconstruction remains available.",
              "- White homoscedastic noise is assumed. Colored noise, drift, unobserved",
              "  time scales, and nearly coincident rates can invalidate count/rate inference.",
              "  Fit reconstruction accuracy alone is not proof of uniquely identified modes.",
              "- Configuration compatibility: `model_selection='bic'` is the new default;",
              "  choose `'aicc'` for AICc. The old `use_aic` flag is retained as the enable",
              "  switch for information-criterion pruning. `disable_als` now skips joint LM",
              "  refinement despite its historical name. Defaults for scales/grid also changed.", "",
              "## All forward measurements", "", *full_table(forward["results"]), "",
              "## All inverse measurements", "", *full_table(inverse["results"]), "",
              "## All component measurements", "",
              "RMSE is against the clean generating signal. Rate error is the largest",
              "relative error after sorting rates; it is omitted when counts disagree.", "",
              *full_table(decomposition["results"], True), "",
              "## Reproduce", "", "```bash", 'python -m pip install ".[benchmark]"',
              "python benchmarks/validate_wavelets.py",
              "python benchmarks/validate_wavelets.py --seed 20261107 --output benchmarks/comparison/validation_fresh.json",
              "python benchmarks/validate_wavelets.py --seed 20261119 --output benchmarks/comparison/validation_holdout.json",
              "python benchmarks/compare.py --section all --repeats 5 --points 64",
              "python benchmarks/compare_extended.py --repeats 5 --points 64",
              "python benchmarks/compare_pairs.py --repeats 21",
              "python benchmarks/render_comparison.py --update-readme", "```", "",
              "Data: [forward](comparison/forward.json), [inverse](comparison/inverse.json),",
              "[inverse high-order follow-up](comparison/inverse_extended.json),",
              "[components](comparison/decomposition.json), [validation](comparison/validation.json),",
              "[fresh validation](comparison/validation_fresh.json),",
              "[additional validation](comparison/validation_holdout.json),",
              "[measurements before the noisy fixes](comparison-before-noisy-fix/validation_fresh.json),",
              "[previous pair-case measurements](comparison-before-pair-optimization/decomposition.json),",
              "[previous algorithm on the fresh seed](comparison-before-pair-optimization/validation_fresh.json).", "",
              "[Longer pair timing check](comparison/pairs.json).", "",
              "## Primary references", "",
              "- [MINPACK QR/Levenberg–Marquardt algorithms](https://www.netlib.org/minpack/)",
              "- [SciPy least_squares: TRF and LM](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)",
              "- [mpmath inversion algorithms and precision](https://mpmath.org/doc/current/calculus/inverselaplace.html)",
              "- [SciPy quad](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html)",
              "- [SciPy quad_vec](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad_vec.html)",
              "- [SciPy impulse](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.impulse.html)", ""]
    (ROOT/"COMPARISON.md").write_text("\n".join(lines), encoding="utf-8")


def update_readme(inverse, decomposition, validation):
    research_section = ROOT / "research" / "README_SECTION.md"
    if research_section.exists():
        path = ROOT.parent / "README.md"
        text = path.read_text(encoding="utf-8")
        start, end = "<!-- comparison:start -->", "<!-- comparison:end -->"
        begin, finish = text.index(start), text.index(end) + len(end)
        replacement = start + "\n" + research_section.read_text(encoding="utf-8").strip() + "\n" + end
        path.write_text(text[:begin] + replacement + text[finish:], encoding="utf-8")
        return
    lines = ["<!-- comparison:start -->", "## Log-wavelet recovery: measured results", "",
             "Persistent log-plane wavelet ridges now initialize a joint rate/amplitude fit.",
             "Noise propagation, BIC selection and refitted deletions control missed and",
             "duplicate modes. The default route does not need the true component count.", "",
             *validation_text(validation), "",
             *fresh_validation_text(), "",
             "The boundary uncertainty check removes an unconstrained extrapolated mode.",
             "Residual MAD calibration reveals a previously suppressed weak component.",
             "Explicit noise settings are respected, and identifiable fast modes remain supported.", "",
             *noisy_fix_summary(), "",
             "![Component speed and weak-mode accuracy](benchmarks/comparison/wavelet_comparison.png)", "",
             *decomposition_summary(decomposition), "",
             "Times are medians of five local runs after warmup. SciPy comparison includes",
             "both TRF and LM automatic-count fits with BIC, plus the original residual-stop",
             "baseline. The table chooses the fastest configuration that recovers the count",
             "and meets the reported reconstruction/rate accuracy criteria. Known-count",
             "SciPy oracle results, every configuration, and individual timings are in the report.", "",
             "### Where the accuracy advantage comes from", "", *weak_summary(decomposition), "",
             "At 1000:1 amplitudes, the residual-stop SciPy baseline can stop at two modes.",
             "Wavelet-guided selection recovers the weak third mode. Stronger SciPy BIC",
             "selection or supplying the count also recovers it, with additional fitting",
             "work. Successful fits generally reach the same float64 accuracy; this is",
             "an automatic model-selection advantage, not universal precision superiority.", "",
             "### Isolating the wavelet contribution", "",
             "The pair-case follow-up caches noise thresholds at crossings, pads filter",
             "boundaries, refines unresolved ridge proposals on both sides, and reuses QR",
             "across damping trials. Count-selection rules and fit tolerances are retained.", "",
             *pair_optimization_summary(decomposition), "",
             *pair_check_summary(), "",
             "The no-wavelet ablation uses rate-grid candidates with the **same new QR/LM",
             "optimizer, noise tests and model selection**. It separates wavelet initialization",
             "from the optimizer rewrite. Full before/after and ablation results are published.", "",
             *ablation_text(decomposition), "",
             "The exact single-mode shortcut is shared by both routes; speed gains on",
             "that case belong to the shortcut rather than to wavelet analysis.", "",
             "![Inverse speed and wavelet ablation](benchmarks/comparison/comparison.png)", "",
             "The ordinary forward quadrature and Weeks inverse paths remain separate from",
             "wavelet decomposition. Their speedups over mpmath are reported independently.",
             "No finite benchmark establishes superiority on all signals: colored noise,",
             "unobserved time scales and nearly coincident modes remain difficult.", "",
             "[Full report, settings, historical measurements, and raw data](benchmarks/COMPARISON.md).", "",
             "```bash", 'python -m pip install ".[benchmark]"',
             "python benchmarks/validate_wavelets.py",
             "python benchmarks/validate_wavelets.py --seed 20261107 --output benchmarks/comparison/validation_fresh.json",
             "python benchmarks/validate_wavelets.py --seed 20261119 --output benchmarks/comparison/validation_holdout.json",
             "python benchmarks/compare.py --section all --repeats 5 --points 64",
             "python benchmarks/compare_extended.py --repeats 5 --points 64",
             "python benchmarks/compare_pairs.py --repeats 21",
             "python benchmarks/render_comparison.py --update-readme", "```", "",
             "<!-- comparison:end -->", ""]
    path = ROOT.parent/"README.md"
    text = path.read_text(encoding="utf-8")
    start, end = "<!-- comparison:start -->", "<!-- comparison:end -->"
    if start in text:
        begin = text.index(start)
        finish = text.index(end, begin)+len(end)
        text = text[:begin]+"\n".join(lines).rstrip()+text[finish:]
    else:
        text = text.replace("## Quick Start", "\n".join(lines)+"\n## Quick Start", 1)
    path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--update-readme", action="store_true")
    args = parser.parse_args()
    forward, inverse, decomposition = [load(section) for section in ("forward", "inverse", "decomposition")]
    extended, validation = load("inverse_extended"), load("validation")
    for dataset in (forward, decomposition, extended):
        for key in ("extension_sha256", "points", "repeats", "script_sha256"):
            if dataset["metadata"][key] != inverse["metadata"][key]:
                raise ValueError(f"Mismatched benchmark metadata: {key}")
    if validation["extension_sha256"] != inverse["metadata"]["extension_sha256"]:
        raise ValueError("Validation used a different extension")
    if load("validation_fresh")["extension_sha256"] != inverse["metadata"]["extension_sha256"]:
        raise ValueError("Fresh validation used a different extension")
    if load("validation_holdout")["extension_sha256"] != inverse["metadata"]["extension_sha256"]:
        raise ValueError("Additional validation used a different extension")
    if load("pairs")["metadata"]["extension_sha256"] != inverse["metadata"]["extension_sha256"]:
        raise ValueError("Pair timing check used a different extension")
    for name, case in extended["results"].items():
        inverse["results"][name]["methods"].update(case["methods"])
    figures(inverse, decomposition)
    wavelet_figure(decomposition)
    write_report(forward, inverse, decomposition, validation)
    if args.update_readme:
        update_readme(inverse, decomposition, validation)
    for path in DATA.glob("*.svg"):
        path.write_text("\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines())+"\n", encoding="utf-8")
    print("Saved comparison report, README section, and PNG/SVG figures")
