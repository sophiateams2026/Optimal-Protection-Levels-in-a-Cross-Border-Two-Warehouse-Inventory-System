"""Reproduce every table, figure and number in Section 6 of the manuscript.

    python reproduce_all.py            # about 30 seconds on a laptop

Outputs (in ./output):
    table1.csv, table2.csv, fig1.png, fig2.png, fig3a.png, fig3b.png, fig4.png,
    results.json             every computed value,
    verification_report.md   reported vs computed values, PASS/FAIL per line.

The exit code is 0 if every check passes and 1 otherwise.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

import hbna  # noqa: E402
import table1_parameters  # noqa: E402
import fig1_accessible_inventory  # noqa: E402
import table2_fig2_expected_cost  # noqa: E402
import fig3_sensitivity  # noqa: E402
import fig4_regimes_exporter  # noqa: E402


def fmt(x: float) -> str:
    """Compact number for the report: 4 decimals, trailing zeros removed."""
    x = 0.0 if abs(x) < 5e-7 else float(x)
    return f"{x:.4f}".rstrip("0").rstrip(".")


def check_internal(spec: dict, results: dict):
    v = results[spec["key"]]
    if "max" in spec:
        return v <= spec["max"], f"<= {spec['max']:g}"
    if "min" in spec:
        return v >= spec["min"], f">= {spec['min']:g}"
    if "abs" in spec:
        return abs(v) <= spec["abs"], f"|x| <= {spec['abs']:g}"
    if "equals_key" in spec:
        w = results[spec["equals_key"]]
        return v == w, f"= {w:g}"
    raise ValueError(spec)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--outdir", default=str(ROOT / "output"))
    args = ap.parse_args(argv)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    steps = [("Table 1", table1_parameters.run),
             ("Fig. 1", fig1_accessible_inventory.run),
             ("Table 2 and Fig. 2", table2_fig2_expected_cost.run),
             ("Fig. 3", fig3_sensitivity.run),
             ("Fig. 4 and exporter scenario", fig4_regimes_exporter.run)]

    results, t_all = {}, time.time()
    for name, fn in steps:
        t0 = time.time()
        r = fn(out)
        r.pop("_table", None)
        results.update(r)
        print(f"  {name:32s} done in {time.time() - t0:5.1f} s")

    spec = json.loads((ROOT / "paper_values.json").read_text(encoding="utf-8"))
    lines, n_fail = [], 0

    lines.append("# Verification report\n")
    lines.append(f"Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC} with Python "
                 f"{platform.python_version()}, NumPy {np.__version__}, "
                 f"Matplotlib {matplotlib.__version__}, hbna {hbna.__version__}.\n")
    lines.append("## A. Numbers reported in the manuscript\n")
    lines.append("| Where | Statement | Reported | Computed | Tolerance | Result |")
    lines.append("|---|---|---:|---:|---:|:---:|")
    for s in spec["values"]:
        got = results[s["key"]]
        ok = abs(got - s["reported"]) <= s["tol"] + 1e-12
        n_fail += not ok
        lines.append(f"| {s['where']} | {s['statement']} | {s['reported']:g} | {fmt(got)} | "
                     f"{s['tol']:g} | {'PASS' if ok else '**FAIL**'} |")

    lines.append("\n## B. Internal consistency checks of the analytical results\n")
    lines.append("| Check | Computed | Criterion | Result |")
    lines.append("|---|---:|---|:---:|")
    for s in spec["internal_checks"]:
        ok, crit = check_internal(s, results)
        n_fail += not ok
        lines.append(f"| {s['statement']} | {results[s['key']]:.3g} | {crit} | "
                     f"{'PASS' if ok else '**FAIL**'} |")

    n_total = len(spec["values"]) + len(spec["internal_checks"])
    summary = (f"{n_total - n_fail} of {n_total} checks passed."
               if n_fail == 0 else f"{n_fail} of {n_total} checks FAILED.")
    lines.insert(2, f"**{summary}**\n")

    (out / "verification_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out / "results.json").write_text(
        json.dumps({k: float(v) for k, v in results.items()}, indent=1, sort_keys=True),
        encoding="utf-8")

    print(f"\n{summary}  ({time.time() - t_all:.1f} s)")
    print(f"Report: {out / 'verification_report.md'}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
