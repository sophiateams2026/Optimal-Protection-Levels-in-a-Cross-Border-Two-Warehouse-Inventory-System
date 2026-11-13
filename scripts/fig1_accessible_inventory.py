"""Fig. 1 (Section 6.2): accessible inventory under DF, P(300) and OF
when a three-week restriction begins at week 2.

The paths come from the time-stepping simulation in hbna/simulation.py,
which applies the issuing rule of Definition 1 directly. The script also
checks the statements made about Fig. 1 in the text of Section 6.2.

Run:  python scripts/fig1_accessible_inventory.py
"""
from __future__ import annotations

import numpy as np

from _common import BLUE, RED, finish, plt
import hbna
from hbna.simulation import simulate

TAU = 2.0
POLICIES = {"DF (R = 0)": 0.0, "P(300)": 300.0, "OF (R = 400)": 400.0}
DT = 0.001


def run(path=None, show: bool = False) -> dict:
    prm = hbna.baseline()
    L, T0 = prm.L, prm.T0
    sims = {name: simulate(prm, R, [TAU], dt=DT, record=True) for name, R in POLICIES.items()}
    t = next(iter(sims.values()))["t"]

    def at(x, time):
        return float(x[int(round(time / DT)), 0])

    res = {}
    keys = {"DF (R = 0)": "DF", "P(300)": "P300", "OF (R = 400)": "OF"}
    during = (t >= TAU) & (t < TAU + L)
    for name, s in sims.items():
        k = keys[name]
        IA, IP = s["I_A"], s["I_D"] + s["I_O"]
        res[f"fig1.{k}.IA_at_onset"] = at(IA, TAU)
        res[f"fig1.{k}.IA_before_restoration"] = at(IA, TAU + L - DT)
        res[f"fig1.{k}.backlog_units"] = float(s["U"][0])
        res[f"fig1.{k}.min_physical_during"] = float(IP[during, 0].min())
        res[f"fig1.{k}.IA_at_restoration"] = at(IA, TAU + L)
        empty = np.where(during & (IA[:, 0] <= 1e-9))[0]
        res[f"fig1.{k}.stockout_week"] = float(t[empty[0]]) if empty.size else float("nan")

    res["fig1.min_physical_all_policies"] = min(
        v for k, v in res.items() if k.endswith("min_physical_during"))

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axvspan(TAU, TAU + L, color="#ececec", zorder=0)
    ax.text(TAU + L / 2, 1040, "OW inaccessible", ha="center")
    common = sims["P(300)"]["I_A"][:, 0]
    before, after = t < TAU, t >= TAU + L
    ax.plot(t[before], common[before], "k", lw=1.4, label="Common path")
    ax.plot(t[after], common[after], "k", lw=1.4)
    styles = {"DF (R = 0)": dict(color="k", ls="--"),
              "P(300)": dict(color=BLUE, ls="-"),
              "OF (R = 400)": dict(color=RED, ls=":", lw=2)}
    for name, s in sims.items():
        ax.plot(t[during], s["I_A"][during, 0], label=name, **styles[name])
    ax.set_xlim(0, T0)
    ax.set_ylim(0, 1100)
    ax.set_xlabel("Time (weeks)")
    ax.set_ylabel(r"Accessible inventory $I_A(t)$")
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.15), frameon=False)
    finish(fig, "fig1.png", path, show)
    return res


if __name__ == "__main__":
    for k, v in run().items():
        print(f"{k:36s} {v:8.2f}")
