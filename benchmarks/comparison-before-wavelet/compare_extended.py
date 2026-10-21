"""Follow-up resolution checks after the initial comparison exposed failures.

Keeps the original matrix intact, adds LaplaceX N512 for every inverse case,
and adds de Hoog degree64 on the delayed step. Run after compare.py, serially.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from compare import (
    HERE,
    Output,
    inverse_cases,
    lx,
    measure,
    metadata,
    mp,
    np,
    reference,
)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--points", type=int, default=64)
    parser.add_argument("--output-dir", type=Path, default=HERE/"comparison")
    args = parser.parse_args()
    if args.repeats < 1 or args.points < 2:
        parser.error("repeats must be positive and points must be at least two")
    results = {}
    t = np.linspace(.1, 5, args.points)
    for name, image, exact, _ in inverse_cases():
        truth = reference(exact, t)
        methods = {"LaplaceX N512": lambda: Output(lx.inverse(image, t, N=512))}
        if name == "delayed_step":
            def dehoog64():
                with mp.workdps(15):
                    values = [float(mp.invertlaplace(image, float(ti),
                                                     method="dehoog", degree=64)) for ti in t]
                return Output(np.asarray(values))
            methods["mpmath deHoog degree64"] = dehoog64
        results[name] = {"times": t.tolist(), "methods": measure(methods, truth, args.repeats)}
        print(f"extended inverse {name} finished", flush=True)
    meta = metadata(args)
    meta["extension_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    meta["purpose"] = "Follow-up higher-resolution settings; initial comparison remains unchanged"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/"inverse_extended.json").write_text(
        json.dumps({"metadata": meta, "results": results}, indent=2, allow_nan=False)+"\n",
        encoding="utf-8")
