"""Time-stepping simulation of the protection-level policy P(R).

This module does not use any formula from the paper. It moves the two
warehouses forward in small time steps and applies the issuing rule of
Definition 1 and the operating assumptions of Assumption 2 literally:

  (i)   issue from the DW while I_D > R;
  (ii)  then from the OW while the OW is accessible and I_O > 0;
  (iii) then from the DW reserve when the OW is empty or inaccessible;
  unmet demand is backlogged and served first when access returns
  (Assumption 2(iii)); stranded OW stock is neither issued nor lost
  (Assumption 2(i)); the next replenishment at T0 clears any backlog.

The simulation is vectorised over many onset times, so averaging over a fine
grid of onsets gives an estimate of E[E_S], E[A_B] and E[U] that can be
compared with Lemma 2 and the closed forms of Lemma 4. The only error
comes from rounding the onset and restoration times to the time grid, so it
is of order dt.
"""
from __future__ import annotations

import numpy as np

from .params import Params


def simulate(prm: Params, R: float, taus, dt: float = 0.001, record: bool = False):
    """Simulate one cycle for each onset in ``taus`` (use np.inf for no restriction).

    Returns a dict with arrays (one entry per onset):
        ES  cumulative stranded inventory  int I_S dt   (unit-weeks)
        AB  cumulative backlog             int B dt     (unit-weeks)
        U   units backordered
        HC  holding cost  int (h_D I_D + h_O I_O) dt over the cycle
    If ``record`` is True, also returns the time grid ``t`` and the paths
    ``I_D``, ``I_O``, ``I_A`` (accessible), ``I_S`` (stranded), ``B``
    with shape (len(t), len(taus)).
    """
    taus = np.atleast_1d(np.asarray(taus, float))
    n = taus.size
    nsteps = int(round(prm.T0 / dt))
    hD = prm.r * (prm.v + prm.dv) + prm.sD
    hO = prm.r * prm.v + prm.sO

    ID = np.full(n, float(prm.QD))
    IO = np.full(n, float(prm.QO))
    B = np.zeros(n)
    ES = np.zeros(n)
    AB = np.zeros(n)
    U = np.zeros(n)
    HC = np.zeros(n)
    was_off = np.zeros(n, bool)

    if record:
        t_rec = np.arange(nsteps + 1) * dt
        rec = {k: np.empty((nsteps + 1, n)) for k in ("I_D", "I_O", "I_A", "I_S", "B")}

    for k in range(nsteps + 1):
        t = k * dt
        off = (taus <= t + 1e-12) & (t + 1e-12 < taus + prm.L)
        on = ~off

        # Access restored: serve the backlog first, from OW then DW.
        restore = on & was_off & (B > 0)
        if restore.any():
            fill = np.where(restore, np.minimum(B, IO), 0.0)
            IO -= fill
            B -= fill
            fill = np.where(restore, np.minimum(B, ID), 0.0)
            ID -= fill
            B -= fill
        was_off = off

        if record:
            rec["I_D"][k] = ID
            rec["I_O"][k] = IO
            rec["I_A"][k] = ID + on * IO
            rec["I_S"][k] = off * IO
            rec["B"][k] = B
        if k == nsteps:
            break

        # Stock levels at the start of the step (for the trapezoid rule below).
        B0 = B.copy()
        H0 = hD * ID + hO * IO
        ES += off * IO * dt          # OW stock is constant while stranded

        need = np.full(n, prm.d * dt)
        take = np.minimum(need, np.maximum(ID - R, 0.0))          # rule (i)
        ID -= take
        need -= take
        take = np.where(on, np.minimum(need, IO), 0.0)            # rule (ii)
        IO -= take
        need -= take
        take = np.minimum(need, ID)                               # rule (iii)
        ID -= take
        need -= take
        B += need                                                 # unmet -> backlog
        U += need
        # Paths are linear within a step, so the trapezoid rule is exact there.
        AB += 0.5 * (B0 + B) * dt
        HC += 0.5 * (H0 + hD * ID + hO * IO) * dt

    # A restriction still running at T0 keeps the remaining OW stock stranded
    # until tau + L (E_S = L * E_0, Eq. (2)); the backlog is cleared at T0.
    tail = np.where(np.isfinite(taus), np.maximum(taus + prm.L - prm.T0, 0.0), 0.0)
    ES += IO * tail

    out = {"ES": ES, "AB": AB, "U": U, "HC": HC}
    if record:
        out.update(rec)
        out["t"] = t_rec
    return out


def onset_grid(prm: Params, n_onsets: int):
    """Midpoint grid of onsets on [0, T0)."""
    return (np.arange(n_onsets) + 0.5) / n_onsets * prm.T0


def expected_costs_sim(prm: Params, R: float, n_onsets: int = 1000, dt: float = 0.001):
    """Average the simulation over a uniform grid of onsets in [0, T0).

    Returns (E[E_S], E[A_B], E[U]), comparable with
    :func:`hbna.model.expected_costs_closed`.
    """
    res = simulate(prm, R, onset_grid(prm, n_onsets), dt)
    return res["ES"].mean(), res["AB"].mean(), res["U"].mean()


def holding_cost_sim(prm: Params, R: float, dt: float = 0.001) -> float:
    """Holding cost of a cycle without disruption, for comparison with Lemma 3."""
    return float(simulate(prm, R, [np.inf], dt)["HC"][0])
