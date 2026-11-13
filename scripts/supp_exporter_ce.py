"""Supplementary analysis (not reported in the manuscript): the exporter
scenario of Section 6.4 with a smaller stranded-inventory cost.

Section 6.4 raises Delta h from 0.12 to 0.13 for the exporter and keeps
c_e = 0.15. Section 3.4 notes that stock stranded at home can often be
redirected, so the write-down part of c_e (delta * v_O) may be smaller for an
exporter. This script re-solves the exporter case with the write-down halved
and removed, and plots R* against c_e for both firm types.

Run:  python scripts/supp_exporter_ce.py
"""
from __future__ import annotations

import csv

import numpy as np

from _common import BLUE, RED, finish, outdir, plt
import hbna
from hbna import model

R_INT = np.arange(0, 401, dtype=float)


def solve(prm) -> dict:
    J = model.J_curve(prm, R_INT)
    k = int(J.argmin())
    c = model.J_components(prm, R_INT[k])
    return {"Rstar": float(R_INT[k]), "J": float(J[k]), "J_DF": float(J[0]), "J_OF": float(J[-1]),
            "premium": c["premium"], "stranded": c["stranded"],
            "backlog": c["backlog_time"] + c["backlog_unit"]}


def run(path=None, show: bool = False) -> dict:
    imp, exp = hbna.baseline(), hbna.exporter()
    scenarios = [
        ("Importer baseline", imp),
        ("Exporter, c_e unchanged (manuscript)", exp),
        ("Exporter, write-down halved", exp.with_(delta=0.0005)),
        ("Exporter, no write-down", exp.with_(delta=0.0)),
    ]
    res, rows = {}, []
    for i, (name, q) in enumerate(scenarios):
        s = solve(q)
        rows.append([name, f"{q.ce:.3f}", f"{q.dh:.3f}", f"{s['Rstar']:.0f}", f"{s['premium']:.2f}",
                     f"{s['stranded']:.2f}", f"{s['backlog']:.2f}", f"{s['J']:.2f}",
                     f"{s['J_DF']:.2f}", f"{s['J_OF']:.2f}",
                     f"{100 * (1 - s['J'] / s['J_DF']):.1f}", f"{100 * (1 - s['J'] / s['J_OF']):.1f}"])
        for k, v in s.items():
            res[f"supp.scenario{i}.{k}"] = v

    header = ["Scenario", "c_e", "Delta h", "R*", "Premium", "Stranded", "Backlog", "J*",
              "J(DF)", "J(OF)", "vs DF %", "vs OF %"]
    with open(outdir(path) / "supp_exporter_ce.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)

    ces = np.linspace(0, 0.30, 301)
    r_imp = [model.Rstar_theorem(imp.with_(ce=c), step=0.05) for c in ces]
    r_exp = [model.Rstar_theorem(exp.with_(ce=c), step=0.05) for c in ces]
    res["supp.exporter.dRstar_dce_0.10_to_0.15"] = (
        solve(exp.with_(ce=0.15))["Rstar"] - solve(exp.with_(ce=0.10))["Rstar"]) / 0.05

    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(ces, r_imp, color=BLUE, lw=1.6, label=r"Importer ($\Delta h$ = 0.12)")
    ax.plot(ces, r_exp, color=RED, lw=1.6, ls="--", label=r"Exporter ($\Delta h$ = 0.13)")
    for _, q in scenarios[1:]:
        ax.plot(q.ce, solve(q)["Rstar"], "s", color=RED, ms=5)
    ax.plot(imp.ce, solve(imp)["Rstar"], "o", color=BLUE, ms=5)
    ax.set_xlabel(r"Stranded-inventory cost $c_e$ (per unit-week)")
    ax.set_ylabel("Optimal reserve R*")
    ax.set_title("p = 0.2, L = 3")
    ax.legend()
    finish(fig, "supp_exporter_ce.png", path, show)

    res["_table"] = (header, rows)
    return res


if __name__ == "__main__":
    out = run()
    header, rows = out.pop("_table")
    widths = [38] + [9] * (len(header) - 1)
    print("".join(f"{h:>{w}s}" if i else f"{h:{w}s}" for i, (h, w) in enumerate(zip(header, widths))))
    for r in rows:
        print("".join(f"{c:>{w}s}" if i else f"{c:{w}s}" for i, (c, w) in enumerate(zip(r, widths))))
    print(f"\nExporter: dR*/dc_e between 0.10 and 0.15 = "
          f"{out['supp.exporter.dRstar_dce_0.10_to_0.15']:.0f} units per 1.0 of c_e")
