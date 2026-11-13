"""Table 2 and Fig. 2 (Section 6.3): expected relevant cost per cycle.

Table 2 is computed three independent ways:
  A. Lemma 2 path by path, averaged over the uniform onset (numerical);
  B. the closed forms of Lemma 4 / Appendix A;
  C. the time-stepping simulation of Definition 1 (hbna/simulation.py).
The script also checks convexity of J(R) (Theorem 1), the right derivative
in Eq. (14) and the minimisers shown in Fig. 2.

Run:  python scripts/table2_fig2_expected_cost.py
"""
from __future__ import annotations

import csv

import numpy as np

from _common import BLUE, RED, finish, outdir, plt
import hbna
from hbna import model
from hbna.simulation import expected_costs_sim

COLUMNS = {"DF": 0.0, "Rstar": 268.0, "OF": 400.0}
ROWS = [("premium", "Resilience premium pi R"),
        ("stranded", "Stranded-inventory cost p ce E[ES]"),
        ("backlog_time", "Backlog-time cost p cb E[AB]"),
        ("backlog_unit", "Per-unit backlog cost p kappa E[U]"),
        ("J", "Expected relevant cost J(R)"),
        ("E_ES", "Memo: E[ES] (unit-weeks)"),
        ("E_U", "Memo: E[U] (units)")]


def run(path=None, show: bool = False, sim_onsets: int = 1000) -> dict:
    prm = hbna.baseline()
    res = {}

    # ---- R* on the integer grid (closed form) ----
    Rgrid = np.arange(0, int(prm.QD) + 1, dtype=float)
    Jv = model.J_curve(prm, Rgrid)
    Rstar = float(Rgrid[Jv.argmin()])
    res["table2.Rstar"] = Rstar
    columns = dict(COLUMNS, Rstar=Rstar)

    # ---- Table 2, three routes ----
    tab = {}
    for col, R in columns.items():
        a = model.J_components(prm, R, "numeric")
        b = model.J_components(prm, R, "closed")
        es, ab, u = expected_costs_sim(prm, R, n_onsets=sim_onsets)
        c = {"premium": prm.pi * R, "stranded": prm.p * prm.ce * es,
             "backlog_time": prm.p * prm.cb * ab, "backlog_unit": prm.p * prm.kappa * u,
             "E_ES": es, "E_U": u}
        c["J"] = c["premium"] + c["stranded"] + c["backlog_time"] + c["backlog_unit"]
        tab[col] = b
        for key, _ in ROWS:
            res[f"table2.{col}.{key}"] = b[key]
            res[f"check.table2.{col}.{key}.numeric_minus_closed"] = a[key] - b[key]
            res[f"check.table2.{col}.{key}.sim_minus_closed"] = c[key] - b[key]

    res["check.table2.max_abs_sim_minus_closed"] = max(
        abs(v) for k, v in res.items() if k.endswith("sim_minus_closed"))
    res["check.table2.max_abs_numeric_minus_closed"] = max(
        abs(v) for k, v in res.items() if k.endswith("numeric_minus_closed"))

    J_DF, J_R, J_OF = tab["DF"]["J"], tab["Rstar"]["J"], tab["OF"]["J"]
    res["text63.saving_vs_DF_pct"] = 100 * (1 - J_R / J_DF)
    res["text63.saving_vs_OF_pct"] = 100 * (1 - J_R / J_OF)
    res["text63.OF_minus_DF"] = J_OF - J_DF
    res["text63.DF_penalty_and_backlog"] = tab["DF"]["backlog_time"] + tab["DF"]["backlog_unit"]

    # ---- Methods A and B agree on the whole grid, for several L ----
    worst = 0.0
    for L in (1.0, 3.0, 4.5, 8.0):
        q = prm.with_(L=L)
        for R in np.linspace(0, prm.QD, 21):
            a = model.expected_costs_numeric(q, R, n=100_000)
            b = model.expected_costs_closed(q, R)
            worst = max(worst, max(abs(x - y) for x, y in zip(a, b)))
    res["check.lemma4.max_abs_diff_numeric_vs_closed"] = worst

    # ---- Theorem 1: convexity and the right derivative, Eq. (14) ----
    p_vals = {0.1: None, 0.2: None, 0.3: None}
    min_d2, max_slope_err = np.inf, 0.0
    for p in p_vals:
        q = prm.with_(p=p)
        J = model.J_curve(q, Rgrid)
        p_vals[p] = J
        min_d2 = min(min_d2, np.diff(J, 2).min())
        h = 1e-4
        for R in np.arange(0.5, prm.QD, 7.0):          # avoid the kink at R = X
            fd = (model.J(q, R + h) - model.J(q, R)) / h
            max_slope_err = max(max_slope_err, abs(fd - float(model.J_right_derivative(q, R))))
        res[f"fig2.p{p}.Rstar"] = float(Rgrid[J.argmin()])
        res[f"fig2.p{p}.min_J_minus_J0"] = float(J.min() - J[0])
        res[f"fig2.p{p}.JOF_minus_J0"] = float(J[-1] - J[0])
    res["check.thm1.min_second_difference"] = float(min_d2)
    res["check.thm1.max_abs_err_Jprime_eq14"] = max_slope_err

    # ---- Lemma 3: budgeted holding cost is linear in R with slope pi ----
    from hbna.simulation import holding_cost_sim
    errs = [holding_cost_sim(prm, R) - model.budgeted_holding_cost(prm, R) for R in (0, 268, 400)]
    res["check.lemma3.max_abs_diff_sim_vs_eq7"] = float(max(abs(e) for e in errs))
    res["check.lemma3.slope_equals_pi"] = (model.budgeted_holding_cost(prm, 400)
                                           - model.budgeted_holding_cost(prm, 0)) / 400 - prm.pi

    # ---- write Table 2 ----
    with open(outdir(path) / "table2.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Cost component", "DF (R=0)", f"P(R*), R*={Rstar:g}", "OF (R=400)"])
        for key, label in ROWS:
            w.writerow([label] + [f"{tab[c][key]:.2f}" for c in columns])

    # ---- Fig. 2 ----
    fig, ax = plt.subplots(figsize=(7, 4.5))
    styles = {0.1: ("k", "--"), 0.2: (BLUE, "-"), 0.3: (RED, ":")}
    for p, (col, ls) in styles.items():
        J = p_vals[p]
        ax.plot(Rgrid, J - J[0], color=col, ls=ls, lw=1.6, label=f"p = {p}")
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_xlabel("Reserve R")
    ax.set_ylabel(r"$J(R) - J(0)$")
    ax.set_xlim(0, 400)
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.13), frameon=False)
    finish(fig, "fig2.png", path, show)
    return res


def print_table(res: dict) -> None:
    print(f"{'Cost component':40s} {'DF':>9s} {'R*':>9s} {'OF':>9s}")
    for key, label in ROWS:
        vals = [res[f"table2.{c}.{key}"] for c in COLUMNS]
        print(f"{label:40s} " + " ".join(f"{v:9.2f}" for v in vals))


if __name__ == "__main__":
    out = run()
    print(f"R* = {out['table2.Rstar']:g}\n")
    print_table(out)
    print()
    sim = max(abs(v) for k, v in out.items() if k.endswith("sim_minus_closed"))
    num = max(abs(v) for k, v in out.items() if k.endswith("numeric_minus_closed"))
    print(f"Table 2: max |Method A - Method B| = {num:.2e}; max |simulation - Method B| = {sim:.3f}")
    for k, v in out.items():
        if not k.startswith(("table2.", "check.table2.")):
            print(f"{k:48s} {v: .6g}")
