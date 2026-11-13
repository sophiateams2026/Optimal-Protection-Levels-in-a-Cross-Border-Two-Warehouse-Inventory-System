"""Fig. 3 (Section 6.4): optimal reserve R* as a function of
(a) the restriction duration L (p = 0.2 and 0.3) and (b) the restriction
probability p (L = 3).

R* is computed from Theorem 1, Eq. (15), and checked against a brute-force
minimisation of J(R) that uses only Lemma 2. The thresholds quoted in
Section 6.4 are solved exactly by bisection. The script also checks
Corollary 2, Corollary 3 (with a random duration), Proposition 2 and
Proposition 3 on parameter grids.

Run:  python scripts/fig3_sensitivity.py
"""
from __future__ import annotations

import numpy as np

from _common import BLUE, RED, bisect, finish, plt
import hbna
from hbna import model


def _quantile_R(prm, Ls, probs, R_grid):
    """Brute-force R* when L is random on a finite support (Corollary 3 check)."""
    J = np.zeros_like(R_grid)
    for L, w in zip(Ls, probs):
        q = prm.with_(L=L)
        J += w * np.array([model.J(q, R) - q.pi * R for R in R_grid])
    J += prm.pi * R_grid
    return float(R_grid[np.isclose(J, J.min(), rtol=0.0, atol=1e-9)].min())


def run(path=None, show: bool = False) -> dict:
    prm = hbna.baseline()
    dhT0 = prm.dh * prm.T0
    res = {}

    # ---- thresholds quoted in Section 6.4 (solved exactly) ----
    q = prm.with_(p=0.2)
    res["text64.p0.2.L_first_positive"] = bisect(
        lambda L: q.p * float(model.M(q.with_(L=L), 0.0)) - dhT0, 1e-6, 4.0)
    res["text64.p0.2.L_first_full"] = bisect(
        lambda L: q.p * float(model.M_left_QD(q.with_(L=L))) - dhT0, 4.0, 8.0)
    res["text64.L3.p_DF_upper"] = dhT0 / float(model.M(prm, 0.0))
    res["text64.L3.p_reaches_300"] = dhT0 / float(model.M(prm, 300 - 1e-9))
    res["text64.L3.Corollary2_max_p_ce_L"] = 1.0 * prm.ce * prm.L      # p <= 1
    res["text64.L3.Corollary2_ds_T0"] = dhT0

    # ---- Theorem 1 vs brute force ----
    worst = 0.0
    for p, L in [(0.2, 3), (0.3, 3), (0.2, 1), (0.2, 4.5), (0.15, 3), (0.1, 3), (0.3, 6)]:
        c = prm.with_(p=p, L=L)
        a, b = model.Rstar_theorem(c), model.Rstar_bruteforce(c)
        res[f"fig3.p{p}.L{L}.Rstar_theorem"] = a
        res[f"fig3.p{p}.L{L}.Rstar_bruteforce"] = b
        worst = max(worst, abs(a - b))
    res["check.thm1.max_abs_diff_theorem_vs_bruteforce"] = worst   # brute-force grid step 0.5

    # ---- curves for Fig. 3 ----
    Ls = np.linspace(0.01, 8, 800)
    curves_L = {p: np.array([model.Rstar_theorem(prm.with_(p=p, L=L), step=0.05) for L in Ls])
                for p in (0.2, 0.3)}
    ps = np.linspace(0, 0.5, 1001)
    curve_p = np.array([model.Rstar_theorem(prm.with_(p=p), step=0.05) for p in ps])
    p_all = np.linspace(0.001, 1.0, 1000)
    res["text64.L3.max_Rstar_any_p_le_1"] = max(
        model.Rstar_theorem(prm.with_(p=p), step=0.05) for p in p_all)

    # ---- Proposition 2: monotonicity of R* ----
    viol = 0
    viol += int((np.diff(curves_L[0.2]) < -1e-9).sum() + (np.diff(curves_L[0.3]) < -1e-9).sum())
    viol += int((np.diff(curve_p) < -1e-9).sum())
    for name, grid, sign in [("ce", np.linspace(0, 1, 101), 1), ("cb", np.linspace(0, 5, 101), 1),
                             ("kappa", np.linspace(0, 15, 101), 1),
                             ("dh", np.linspace(0.01, 0.4, 101), -1)]:
        r = np.array([model.Rstar_theorem(prm.with_(**{name: x}), step=0.05) for x in grid])
        viol += int((sign * np.diff(r) < -1e-9).sum())
    res["check.prop2.monotonicity_violations"] = viol

    # ---- Proposition 3: V(p) = J(0) - J(R*) is 0 below the threshold, nondecreasing, convex ----
    Rg = np.arange(0, 401, dtype=float)
    pv = np.linspace(0, 0.5, 201)
    V = []
    for p in pv:
        J = model.J_curve(prm.with_(p=p), Rg)
        V.append(J[0] - J.min())
    V = np.array(V)
    thr = res["text64.L3.p_DF_upper"]
    res["check.prop3.max_V_below_threshold"] = float(V[pv <= thr].max())
    res["check.prop3.min_first_difference"] = float(np.diff(V).min())
    res["check.prop3.min_second_difference"] = float(np.diff(V, 2).min())

    # ---- Corollary 3 with a random duration (cb = 0, L uniform on {1,...,5}) ----
    c3 = prm.with_(cb=0.0, p=0.3)
    support, probs = np.arange(1, 6, dtype=float), np.full(5, 0.2)
    theta = (c3.dh * c3.T0 / c3.p - c3.ce * float(support @ probs)) / c3.kappa
    cdf = np.cumsum(probs)
    l_q = support[np.argmax(cdf >= 1 - theta - 1e-12)]
    res["check.cor3.theta"] = theta
    res["check.cor3.Rstar_formula"] = min(c3.d * l_q, c3.QD)
    res["check.cor3.Rstar_bruteforce"] = _quantile_R(c3, support, probs, Rg)

    # ---- Fig. 3(a) ----
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for p, (col, ls) in {0.2: (BLUE, "-"), 0.3: (RED, ":")}.items():
        ax.plot(Ls, curves_L[p], color=col, ls=ls, lw=1.8, label=f"p = {p}")
    ax.set_xlabel("Restriction duration L (weeks)")
    ax.set_ylabel("Optimal reserve R*")
    ax.set_xlim(0, 8)
    ax.set_ylim(-5, 450)
    ax.legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.13), frameon=False)
    finish(fig, "fig3a.png", path, show)

    # ---- Fig. 3(b) ----
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.plot(ps, curve_p, color=BLUE, lw=1.8, label="L = 3")
    ax.set_xlabel("Restriction probability p")
    ax.set_ylabel("Optimal reserve R*")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(-5, 450)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.13), frameon=False)
    finish(fig, "fig3b.png", path, show)
    return res


if __name__ == "__main__":
    for k, v in run().items():
        print(f"{k:50s} {v: .6g}")
