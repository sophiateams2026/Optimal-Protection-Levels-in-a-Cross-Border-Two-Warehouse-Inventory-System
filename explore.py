"""Try your own parameter values.

Change any input of Table 1 and the script re-solves the model, prints a
Table 2 style cost breakdown, and cross-checks the answer three ways:
  A. Lemma 2 path by path, averaged over the onset (numerical);
  B. the closed forms of Lemma 4 / Appendix A;
  C. a time-stepping simulation of the issuing rule in Definition 1.
R* from Theorem 1 is also compared with a brute-force search over R.

Examples
--------
  python explore.py                              # Table 1 baseline
  python explore.py --list                       # all parameter names
  python explore.py --set L=4 lam=0.03           # longer, more frequent restrictions
  python explore.py --set p=0.3                  # set p directly instead of lam*T0
  python explore.py --set dv=5                   # exporter case of Section 6.4
  python explore.py --set ce=0.10 cb=2.5 kappa=8 # set cost rates directly
  python explore.py --R 150 300                  # also evaluate reserves 150 and 300
  python explore.py --config my_params.json      # read values from a JSON file
  python explore.py --sweep L 0.5 8 31           # R* and costs over a grid of L
  python explore.py --sweep kappa 0 20 41 --set L=2

Sweep tables and figures go to output/explore/.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402

from _common import BLUE, RED, finish, plt  # noqa: E402
import hbna  # noqa: E402
from hbna import model  # noqa: E402
from hbna.simulation import expected_costs_sim  # noqa: E402

OUT = ROOT / "output" / "explore"

DESCRIPTIONS = {
    "Q": "replenishment quantity (units)",
    "d": "demand rate (units/week)",
    "QD": "allocation to demand-side warehouse, 0 < QD < Q",
    "v": "landed unit cost (v_O in the export case)",
    "dv": "carrying-amount gap v_D - v_O (0 = importer, 5 = exporter in Sec. 6.4)",
    "P": "selling price",
    "r": "cost of capital per week",
    "sD": "storage rate, demand-side warehouse (per unit-week)",
    "sO": "storage rate, origin-side warehouse (per unit-week)",
    "compliance": "compliance, insurance and bonded fees in c_e",
    "delta": "write-down hazard rate in c_e (c_e = compliance + delta*v)",
    "penalty": "contractual penalty in c_b (c_b = r v + r P + penalty)",
    "kappa": "per-unit backlog cost",
    "lam": "restriction arrival rate per week (p = lam * T0)",
    "L": "restriction duration in weeks",
    "p": "[shortcut] restriction probability per cycle, overrides lam",
    "ce": "[shortcut] stranded-inventory cost, overrides compliance/delta",
    "cb": "[shortcut] backlog-time cost, overrides r v + r P + penalty",
    "dh": "[shortcut] holding-cost differential, overrides sD - sO + r dv",
}


# ---------------------------------------------------------------------------
# Building and checking a parameter set
# ---------------------------------------------------------------------------
def make_params(base: hbna.Params | None = None, **values) -> hbna.Params:
    """Params from keyword values; raises ValueError for unknown names or invalid values."""
    base = base or hbna.baseline()
    fields = {f.name for f in dataclasses.fields(hbna.Params)}
    unknown = [k for k in values if k not in fields and k not in ("p", "ce", "cb", "dh")]
    if unknown:
        raise ValueError(f"unknown parameter(s): {unknown}. Use --list to see the names.")
    prm = base.with_(**{k: float(v) for k, v in values.items()})
    validate(prm)
    return prm


def validate(prm: hbna.Params) -> list[str]:
    """Raise on values outside the model; return warnings for values that strain assumptions."""
    errors = []
    if prm.Q <= 0 or prm.d <= 0:
        errors.append("Q and d must be positive")
    if not 0 < prm.QD < prm.Q:
        errors.append("need 0 < QD < Q (both warehouses hold stock)")
    if prm.L <= 0:
        errors.append("L must be positive")
    if not 0 <= prm.p <= 1:
        errors.append(f"p = {prm.p:.3f} must lie in [0, 1]")
    for name in ("ce", "cb", "kappa"):
        if getattr(prm, name) < 0:
            errors.append(f"{name} must be nonnegative")
    if errors:
        raise ValueError("; ".join(errors))
    warnings = []
    if prm.p > 0.3:
        warnings.append(f"p = {prm.p:.2f} is large. Assumption 1 (at most one restriction per "
                        "cycle) becomes a rough approximation.")
    if prm.L > prm.T0:
        warnings.append(f"L = {prm.L:g} exceeds the cycle length T0 = {prm.T0:g}. The model "
                        "clears backlog at T0 and ignores carry-over (Assumption 2(iii)).")
    if prm.dh <= 0:
        warnings.append("Delta h <= 0: origin-side storage is not cheaper, so OF is optimal "
                        "for any p > 0 (Corollary 2).")
    return warnings


def load_config(path: str) -> dict:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return {k: v for k, v in data.items() if not k.startswith("_")}


def parse_set(items) -> dict:
    out = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"--set expects name=value, got '{item}'")
        k, v = item.split("=", 1)
        out[k.strip()] = float(v)
    return out


# ---------------------------------------------------------------------------
# Analysis of one parameter set
# ---------------------------------------------------------------------------
def pct_below(J: float, ref: float) -> float:
    """Percentage by which J lies below ref; NaN when ref <= 0.

    J(R) is defined net of terms that do not depend on R (Eq. (9)), so it can be
    zero or negative (for example when p = 0 or Delta h < 0). A percentage is
    then not meaningful.
    """
    return 100 * (1 - J / ref) if ref > 1e-12 else float("nan")


def solve(prm: hbna.Params) -> dict:
    """R* (Theorem 1 and integer grid) and the costs of DF, R*, OF."""
    R_cont = model.Rstar_theorem(prm, step=prm.QD / 40000)
    grid = np.arange(0, int(np.floor(prm.QD)) + 1, dtype=float)
    if grid[-1] != prm.QD:
        grid = np.append(grid, prm.QD)
    J = model.J_curve(prm, grid)
    k = int(J.argmin())
    return {"Rstar_theorem": R_cont, "Rstar_grid": float(grid[k]), "Jstar": float(J[k]),
            "J_DF": float(J[0]), "J_OF": float(J[-1]), "regime": model.regime(prm),
            "grid": grid, "J": J}


def cross_check(prm: hbna.Params, reserves, n_onsets: int = 500, steps: int = 5000):
    """Expected costs of each reserve by three independent routes.

    A  Lemma 2 averaged over 100,000 onsets (numerical integration)
    B  closed forms of Lemma 4 / Appendix A
    C  time-stepping simulation of Definition 1 on ``n_onsets`` onsets

    A vs B tests the closed forms. C is compared with Lemma 2 on the *same*
    onsets (A'), so the comparison tests the inventory dynamics and is not
    blurred by quadrature error from using fewer onsets.
    Gaps are reported relative to the size of the cost terms of J(R).
    """
    from hbna.simulation import onset_grid, simulate

    taus = onset_grid(prm, n_onsets)
    w = np.array([prm.p * prm.ce, prm.p * prm.cb, prm.p * prm.kappa])   # cost weights
    rows = []
    for R in reserves:
        a = np.array(model.expected_costs_numeric(prm, R, n=100_000))
        b = np.array(model.expected_costs_closed(prm, R))
        a2 = np.array(model.expected_costs_on_onsets(prm, R, taus))
        sim = simulate(prm, R, taus, dt=prm.T0 / steps)
        c = np.array([sim["ES"].mean(), sim["AB"].mean(), sim["U"].mean()])
        scale = max(abs(prm.pi * R) + w @ np.abs(b), 1.0)
        rows.append({"R": R, "closed": b,
                     "A_vs_B": float(abs(w @ (a - b)) / scale),
                     "C_vs_A": float(abs(w @ (c - a2)) / scale)})
    return rows


def analyse(prm: hbna.Params, extra_R=(), check: bool = True,
            show: bool = False, quiet: bool = False) -> dict:
    """Solve the model for ``prm``, print a report and return the numbers."""
    warnings = validate(prm)
    s = solve(prm)
    columns = [("DF", 0.0), ("R*", s["Rstar_grid"]), ("OF", prm.QD)]
    columns += [(f"R={R:g}", float(R)) for R in extra_R]
    comps = {name: model.J_components(prm, R) for name, R in columns}
    bf = model.Rstar_bruteforce(prm, np.linspace(0, prm.QD, 801), n=4000)

    out = {"params": prm, "solution": s, "components": comps, "Rstar_bruteforce": bf,
           "warnings": warnings}

    if not quiet:
        p = prm
        print("Parameters")
        print(f"  Q={p.Q:g}  d={p.d:g}  QD={p.QD:g}  QO={p.QO:g}  T0={p.T0:g}  L={p.L:g}  X=dL={p.X:g}")
        print(f"  Delta s={p.ds:.4g}  Delta v={p.dv:g}  Delta h={p.dh:.4g}  pi=Delta h*QO/d={p.pi:.4g}")
        print(f"  c_e={p.ce:.4g}  c_b={p.cb:.4g}  kappa={p.kappa:g}  p={p.p:.4g}")
        for w in warnings:
            print(f"  NOTE: {w}")
        print("\nOptimal protection level")
        print(f"  Regime (Theorem 1)            : {s['regime']}")
        print(f"  R* from Eq. (15)              : {s['Rstar_theorem']:.2f}")
        print(f"  R* on integer grid            : {s['Rstar_grid']:g}")
        print(f"  R* by brute force (Lemma 2)   : {bf:.1f}   (grid step {prm.QD / 800:.2g})")
        print(f"  J(R*) = {s['Jstar']:.2f};  J(DF) = {s['J_DF']:.2f};  J(OF) = {s['J_OF']:.2f}")
        a, b = pct_below(s["Jstar"], s["J_DF"]), pct_below(s["Jstar"], s["J_OF"])
        if np.isfinite(a) and np.isfinite(b):
            print(f"  J(R*) is {a:.1f}% below DF and {b:.1f}% below OF")

        print("\nExpected relevant cost per cycle (Table 2 layout)")
        labels = [("premium", "Resilience premium pi R"),
                  ("stranded", "Stranded-inventory cost"),
                  ("backlog_time", "Backlog-time cost"),
                  ("backlog_unit", "Per-unit backlog cost"),
                  ("J", "Expected relevant cost J(R)"),
                  ("E_ES", "Memo: E[ES] (unit-weeks)"),
                  ("E_U", "Memo: E[U] (units)")]
        head = "".join(f"{f'{n} ({R:g})' if n in ('R*',) else n:>13s}" for n, R in columns)
        print(f"  {'':30s}{head}")
        for key, lab in labels:
            print(f"  {lab:30s}" + "".join(f"{comps[n][key]:13.2f}" for n, _ in columns))

    if check:
        rows = cross_check(prm, [R for _, R in columns])
        out["cross_check"] = rows
        ok = all(r["A_vs_B"] < 1e-6 and r["C_vs_A"] < 1e-3 for r in rows) and \
            abs(s["Rstar_theorem"] - bf) <= prm.QD / 800 + 1e-9
        out["all_checks_pass"] = ok
        if not quiet:
            print("\nCross-check: relative gap in J(R) between computation routes")
            print(f"  {'Reserve':>10s} {'A vs B':>12s} {'C vs A':>12s}")
            for r in rows:
                print(f"  {r['R']:10.2f} {r['A_vs_B']:12.2e} {r['C_vs_A']:12.2e}")
            print("  A = Lemma 2 averaged over the onset, B = closed forms of Lemma 4,")
            print("  C = time-stepping simulation of Definition 1 (compared with Lemma 2 on the same onsets)")
            print(f"  Result: {'PASS' if ok else 'CHECK FAILED'} "
                  "(A vs B < 1e-6; C vs A < 0.1%, the simulation error is of order dt; "
                  "R* from Eq. (15) within one brute-force grid step)")

    return out


# ---------------------------------------------------------------------------
# Sweep one parameter
# ---------------------------------------------------------------------------
def sweep(prm: hbna.Params, name: str, lo: float, hi: float, n: int,
          show: bool = False, quiet: bool = False):
    xs = np.linspace(lo, hi, int(n))
    rows = []
    for x in xs:
        q = make_params(prm, **{name: x})
        s = solve(q)
        rows.append([x, s["regime"], s["Rstar_theorem"], s["Rstar_grid"], s["Jstar"],
                     s["J_DF"], s["J_OF"], pct_below(s["Jstar"], s["J_DF"]),
                     pct_below(s["Jstar"], s["J_OF"])])
    header = [name, "regime", "R*_Eq15", "R*_grid", "J*", "J(DF)", "J(OF)", "vs_DF_%", "vs_OF_%"]
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / f"sweep_{name}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    if not quiet:
        print(f"{name:>9s} {'regime':>8s} {'R*':>9s} {'J*':>9s} {'J(DF)':>9s} {'J(OF)':>9s} "
              f"{'vs DF%':>7s} {'vs OF%':>7s}")
        for r in rows:
            print(f"{r[0]:9.4g} {r[1]:>8s} {r[2]:9.2f} {r[4]:9.2f} {r[5]:9.2f} {r[6]:9.2f} "
                  f"{r[7]:7.1f} {r[8]:7.1f}")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(xs, [r[2] for r in rows], color=BLUE, lw=1.6)
    ax[0].set_xlabel(name)
    ax[0].set_ylabel("Optimal reserve R*")
    ax[1].plot(xs, [r[4] for r in rows], color=BLUE, lw=1.6, label="J(R*)")
    ax[1].plot(xs, [r[5] for r in rows], "k--", lw=1.2, label="J(DF)")
    ax[1].plot(xs, [r[6] for r in rows], color=RED, ls=":", lw=1.6, label="J(OF)")
    ax[1].set_xlabel(name)
    ax[1].set_ylabel("Expected relevant cost")
    ax[1].legend()
    f = finish(fig, f"sweep_{name}.png", OUT, show)
    if not quiet:
        print(f"\nTable: {OUT / f'sweep_{name}.csv'}\nFigure: {f}")
    return header, rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--set", nargs="+", metavar="NAME=VALUE", help="change parameter values")
    ap.add_argument("--config", help="JSON file with parameter values (see params_baseline.json)")
    ap.add_argument("--R", nargs="+", type=float, default=[], help="extra reserves to evaluate")
    ap.add_argument("--sweep", nargs=4, metavar=("NAME", "FROM", "TO", "N"),
                    help="solve over N values of one parameter")
    ap.add_argument("--no-check", action="store_true", help="skip the three-way cross-check")
    ap.add_argument("--list", action="store_true", help="list parameter names and baseline values")
    args = ap.parse_args(argv)

    if args.list:
        base = hbna.baseline()
        for k, desc in DESCRIPTIONS.items():
            val = getattr(base, k)
            print(f"  {k:11s} {val:<8g} {desc}")
        return 0

    values = load_config(args.config) if args.config else {}
    values.update(parse_set(args.set))
    try:
        prm = make_params(**values)
    except ValueError as e:
        print(f"Error: {e}")
        return 2

    if args.sweep:
        name, lo, hi, n = args.sweep
        try:
            sweep(prm, name, float(lo), float(hi), int(n))
        except ValueError as e:
            print(f"Error: {e}")
            return 2
        return 0

    for R in args.R:
        if not 0 <= R <= prm.QD:
            print(f"Error: reserve R = {R:g} must lie in [0, QD = {prm.QD:g}]")
            return 2
    res = analyse(prm, extra_R=args.R, check=not args.no_check)
    return 0 if res.get("all_checks_pass", True) else 1


if __name__ == "__main__":
    sys.exit(main())
