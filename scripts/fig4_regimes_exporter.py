"""Fig. 4 (Section 6.4): policy regimes of Theorem 1 at p = 0.2, and the
exporter scenario (v_O = 50, v_D = 55, so Delta h = 0.13).

Upper boundary  Delta h = p M(0)/T0      : DF optimal above it.
Lower boundary  Delta h = p M(Q_D^-)/T0  : OF optimal below it.
The boundaries are compared with a brute-force regime map built from Lemma 2
alone, point by point on a (L, Delta h) grid.

Run:  python scripts/fig4_regimes_exporter.py
"""
from __future__ import annotations

import numpy as np

from _common import BLUE, RED, bisect, finish, plt
import hbna
from hbna import model


def run(path=None, show: bool = False) -> dict:
    prm = hbna.baseline().with_(p=0.2)
    res = {}

    # ---- boundary values at selected durations ----
    for L in (0.25, 1, 3, 3.999, 4.0, 6, 8):
        up, lo = model.regime_boundaries(prm, [L])
        res[f"fig4.L{L}.upper"] = float(up[0])
        res[f"fig4.L{L}.lower"] = float(lo[0])

    # ---- partial-protection range at the baseline differential 0.12 ----
    res["text64.partial_region_L_from"] = bisect(
        lambda L: model.regime_boundaries(prm, [L])[0][0] - 0.12, 1e-6, 4.0)
    res["text64.partial_region_L_to"] = bisect(
        lambda L: model.regime_boundaries(prm, [L])[1][0] - 0.12, 4.0, 8.0)
    res["text64.lower_boundary_jump_at_L4"] = float(
        model.regime_boundaries(prm, [4.0])[1][0] - model.regime_boundaries(prm, [4.0 - 1e-9])[1][0])

    # ---- brute-force regime map (Lemma 2 only) vs Theorem 1 ----
    Rg = np.linspace(0, prm.QD, 401)
    mismatch = total = 0
    for L in np.linspace(0.3, 8, 40):
        up, lo = (float(x[0]) for x in model.regime_boundaries(prm, [L]))
        for dh in np.linspace(0.005, 0.25, 40):
            if min(abs(dh - up), abs(dh - lo)) < 0.002:     # skip points on a boundary
                continue
            rs = model.Rstar_bruteforce(prm.with_(L=L, dh=dh), Rg, n=2000)
            brute = "DF" if rs == 0 else ("OF" if rs == prm.QD else "partial")
            theory = "DF" if dh >= up else ("OF" if dh <= lo else "partial")
            total += 1
            mismatch += brute != theory
    res["check.fig4.regime_map_points"] = total
    res["check.fig4.regime_map_mismatches"] = mismatch

    # ---- exporter scenario vs importer baseline (Section 6.4) ----
    Rint = np.arange(0, 401, dtype=float)
    for name, q in (("importer", hbna.baseline()), ("exporter", hbna.exporter())):
        J = model.J_curve(q, Rint)
        k = int(J.argmin())
        res[f"exporter.{name}.dh"] = q.dh
        res[f"exporter.{name}.Rstar"] = float(Rint[k])
        res[f"exporter.{name}.Jstar"] = float(J[k])
        res[f"exporter.{name}.J_DF"] = float(J[0])
        res[f"exporter.{name}.J_OF"] = float(J[-1])
        res[f"exporter.{name}.saving_vs_DF_pct"] = 100 * (1 - J[k] / J[0])
        res[f"exporter.{name}.saving_vs_OF_pct"] = 100 * (1 - J[k] / J[-1])
        res[f"exporter.{name}.regime"] = {"DF": 0, "partial": 1, "OF": 2}[model.regime(q)]

    # ---- Fig. 4 ----
    Lu = np.linspace(0.25, 8, 800)
    L1 = np.linspace(0, 4, 200, endpoint=False)
    L2 = np.linspace(4, 8, 400)
    up_u, _ = model.regime_boundaries(prm, Lu)
    _, lo1 = model.regime_boundaries(prm, L1)
    _, lo2 = model.regime_boundaries(prm, L2)
    dh_exp = hbna.exporter().dh
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.plot(Lu, up_u, "k-", lw=1.5, label=r"DF boundary $pM(0)/T_0$")
    ax.plot(L1, lo1, "k--", lw=1.5, label=r"OF boundary $pM(Q_D^-)/T_0$")
    ax.plot(L2, lo2, "k--", lw=1.5)
    ax.axhline(0.12, color="grey", lw=0.6, ls=":")
    ax.plot(3, 0.12, "o", color=BLUE, label=r"Importer baseline ($\Delta v=0,\ \Delta h=0.12$)")
    ax.plot(3, dh_exp, "s", color=RED, label=r"Exporter ($\Delta v=5,\ \Delta h=0.13$)")
    ax.text(1, 0.24, r"DF optimal ($R^*=0$)")
    ax.text(1.3, 0.07, r"Partial reserve ($0<R^*<Q_D$)")
    ax.text(5.6, 0.06, r"OF optimal ($R^*=Q_D$)")
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 0.26)
    ax.set_xlabel("Restriction duration L (weeks)")
    ax.set_ylabel("Holding-cost differential $\\Delta h$\n(= $\\Delta s$ when $\\Delta v$ = 0)")
    ax.legend(fontsize=8, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.2), frameon=False)
    finish(fig, "fig4.png", path, show)
    return res


if __name__ == "__main__":
    for k, v in run().items():
        print(f"{k:44s} {v: .6g}")
