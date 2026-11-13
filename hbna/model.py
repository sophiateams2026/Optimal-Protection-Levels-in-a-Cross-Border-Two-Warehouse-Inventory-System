"""Analytical results of the manuscript.

Equation and result numbers refer to the manuscript:

    Eq. (5)-(6)   Lemma 2, state at onset and path-wise backlog
    Eq. (7)-(8)   Lemma 3, budgeted holding cost and resilience premium
    Eq. (9)       expected relevant cost J(R)
    Eq. (10)      psi(c)
    Eq. (11)-(12) Lemma 4, expected exposure and backlog
                  (closed forms of E[U|L] and E[A_B|L] in Appendix A)
    Eq. (13)      marginal avoided loss M(R)
    Eq. (14)-(15) Theorem 1, optimal protection level

Two routes to the expected costs are provided so that each can be checked
against the other:

* ``expected_costs_numeric``: Lemma 2 applied path by path and averaged over
  the uniform onset tau in [0, T0) by midpoint quadrature (Method A).
* ``expected_costs_closed``: the closed forms of Lemma 4 / Appendix A (Method B).

The restriction duration L is deterministic, as in the baseline of Section 6.
All functions take a :class:`hbna.params.Params` object as first argument.
In the export case Delta h replaces Delta s throughout (Section 3.4), which is
what ``prm.dh`` and ``prm.pi`` return.
"""
from __future__ import annotations

import numpy as np

from .params import Params


# ---------------------------------------------------------------------------
# Lemma 2: state at onset and path-wise outcomes for a given onset tau
# ---------------------------------------------------------------------------
def onset_state(prm: Params, R, tau):
    """Eq. (5): demand-side stock I_D(tau-) and stranded stock E_0 at the onset."""
    tau = np.asarray(tau, float)
    left = prm.Q - prm.d * tau
    ID = np.minimum(np.maximum(prm.QD - prm.d * tau, R), left)
    E0 = left - ID
    return ID, E0


def pathwise_outcome(prm: Params, R, tau):
    """Eq. (6): backordered units U, cumulative backlog A_B and exposure E_S."""
    ID, E0 = onset_state(prm, R, tau)
    left = prm.Q - prm.d * np.asarray(tau, float)
    U = np.maximum(np.minimum(prm.X, left) - ID, 0.0)
    AB = U ** 2 / (2 * prm.d)
    ES = prm.L * E0
    return U, AB, ES


# ---------------------------------------------------------------------------
# Expected costs over the uniform onset (Assumption 1)
# ---------------------------------------------------------------------------
def expected_costs_numeric(prm: Params, R: float, n: int = 200_000):
    """Method A: average Lemma 2 over tau ~ U[0, T0) with n midpoints.

    Returns (E[E_S], E[A_B], E[U]).
    """
    tau = (np.arange(n) + 0.5) / n * prm.T0
    U, AB, ES = pathwise_outcome(prm, R, tau)
    return ES.mean(), AB.mean(), U.mean()


def expected_costs_on_onsets(prm: Params, R: float, taus):
    """Lemma 2 averaged over a given set of onsets. Returns (E[E_S], E[A_B], E[U])."""
    U, AB, ES = pathwise_outcome(prm, R, np.asarray(taus, float))
    return ES.mean(), AB.mean(), U.mean()


def _Phi(prm: Params, k: int, c: float) -> float:
    """Phi_k(c) = int_0^c (min{u, Q_O})^k du (Appendix A)."""
    QO = prm.QO
    if c <= 0:
        return 0.0
    if c <= QO:
        return c ** (k + 1) / (k + 1)
    return QO ** (k + 1) / (k + 1) + QO ** k * (c - QO)


def _Psi(prm: Params, k: int, c: float) -> float:
    """Psi_k(c) = int_0^{Q_O} (min{c, z})^k dz (Appendix A)."""
    QO = prm.QO
    if c <= 0:
        return 0.0
    if c >= QO:
        return QO ** (k + 1) / (k + 1)
    return c ** (k + 1) / (k + 1) + c ** k * (QO - c)


def expected_E0(prm: Params, R: float) -> float:
    """Eq. (11): E[E_0] = (Q_O/Q)(Q_D - R + Q_O/2)."""
    return prm.QO / prm.Q * (prm.QD - R + prm.QO / 2)


def expected_costs_closed(prm: Params, R: float):
    """Method B: closed forms of Lemma 4 (Appendix A).

    Returns (E[E_S], E[A_B], E[U]).
    """
    X, Q, d, QD = prm.X, prm.Q, prm.d, prm.QD
    c = X - R
    EU = (_Phi(prm, 1, c) + _Psi(prm, 1, c) - _Phi(prm, 1, X - QD)) / Q
    EAB = (_Phi(prm, 2, c) + _Psi(prm, 2, c) - _Phi(prm, 2, X - QD)) / (2 * d * Q)
    EES = prm.L * expected_E0(prm, R)
    return EES, EAB, EU


def J_components(prm: Params, R: float, method: str = "closed", n: int = 200_000) -> dict:
    """Eq. (9): expected relevant cost of P(R), split as in Table 2."""
    if method == "closed":
        EES, EAB, EU = expected_costs_closed(prm, R)
    elif method == "numeric":
        EES, EAB, EU = expected_costs_numeric(prm, R, n)
    else:
        raise ValueError("method must be 'closed' or 'numeric'")
    out = {
        "premium": prm.pi * R,
        "stranded": prm.p * prm.ce * EES,
        "backlog_time": prm.p * prm.cb * EAB,
        "backlog_unit": prm.p * prm.kappa * EU,
        "E_ES": EES,
        "E_AB": EAB,
        "E_U": EU,
    }
    out["J"] = out["premium"] + out["stranded"] + out["backlog_time"] + out["backlog_unit"]
    return out


def J(prm: Params, R: float, method: str = "closed") -> float:
    """Expected relevant cost J(R), Eq. (9)."""
    return J_components(prm, R, method)["J"]


def J_curve(prm: Params, R_grid, method: str = "closed") -> np.ndarray:
    return np.array([J(prm, R, method) for R in R_grid])


# ---------------------------------------------------------------------------
# Lemma 3: budgeted holding cost of a cycle without disruption
# ---------------------------------------------------------------------------
def budgeted_holding_cost(prm: Params, R: float) -> float:
    """Eq. (7), with s_k replaced by h_k = r v_k + s_k (Section 3.4).

    For an importer (Delta v = 0) this is exactly Eq. (7):
    r v Q^2/(2d) + [s_D Q_D^2 + s_O (Q_O^2 + 2 Q_O Q_D)]/(2d) + (Delta s Q_O/d) R.
    """
    Q, d, QD, QO = prm.Q, prm.d, prm.QD, prm.QO
    hO = prm.r * prm.v + prm.sO
    hD = prm.r * (prm.v + prm.dv) + prm.sD
    return (hD * QD ** 2 + hO * (QO ** 2 + 2 * QO * QD)) / (2 * d) + prm.dh * QO / d * R


# ---------------------------------------------------------------------------
# Eq. (10), Eq. (13) and Theorem 1
# ---------------------------------------------------------------------------
def psi(prm: Params, c):
    """Eq. (10)."""
    c = np.asarray(c, float)
    QO = prm.QO
    return np.where(c < 0, 0.0, np.where(c <= QO, c - c ** 2 / (2 * QO), QO / 2))


def M(prm: Params, R):
    """Eq. (13) for a deterministic duration L: marginal avoided loss of the reserve."""
    R = np.asarray(R, float)
    X = prm.X
    return prm.ce * prm.L + prm.kappa * (X > R) + prm.cb / prm.d * psi(prm, X - R)


def M_left_QD(prm: Params):
    """M(Q_D^-): left limit at Q_D, where P(X > R) tends to 1{X >= Q_D}."""
    X = prm.X
    return prm.ce * prm.L + prm.kappa * (X >= prm.QD) + prm.cb / prm.d * psi(prm, X - prm.QD)


def J_right_derivative(prm: Params, R):
    """Eq. (14): J'(R) = (Q_O/Q)(Delta h T0 - p M(R))."""
    return prm.QO / prm.Q * (prm.dh * prm.T0 - prm.p * M(prm, R))


def Rstar_theorem(prm: Params, step: float = 0.01) -> float:
    """Eq. (15): R* = min{R in [0, Q_D] : p M(R) <= Delta h T0}, else Q_D."""
    grid = np.arange(0.0, prm.QD + step / 2, step)
    ok = prm.p * M(prm, grid) <= prm.dh * prm.T0
    return float(grid[ok.argmax()]) if ok.any() else float(prm.QD)


def Rstar_bruteforce(prm: Params, R_grid=None, n: int = 4000) -> float:
    """Smallest minimiser of J(R) on a grid, using Lemma 2 path by path.

    This route does not use Theorem 1 or the closed forms.
    """
    if R_grid is None:
        R_grid = np.linspace(0, prm.QD, 801)
    R_grid = np.asarray(R_grid, float)
    tau = (np.arange(n) + 0.5) / n * prm.T0
    U, AB, ES = pathwise_outcome(prm, R_grid[:, None], tau[None, :])
    Jv = prm.pi * R_grid + prm.p * (prm.ce * ES.mean(1) + prm.cb * AB.mean(1)
                                    + prm.kappa * U.mean(1))
    return float(R_grid[np.isclose(Jv, Jv.min(), rtol=0.0, atol=1e-9)].min())


def regime(prm: Params) -> str:
    """Policy regime of Theorem 1: 'DF', 'OF' or 'partial'."""
    if prm.p * M(prm, 0.0) <= prm.dh * prm.T0:
        return "DF"
    if prm.p * M_left_QD(prm) > prm.dh * prm.T0:
        return "OF"
    return "partial"


def regime_boundaries(prm: Params, L):
    """Fig. 4 boundaries in the (L, Delta h) plane at the given p.

    Returns (upper, lower): DF is optimal for Delta h >= upper,
    OF for Delta h <= lower, a partial reserve in between.
    """
    L = np.atleast_1d(np.asarray(L, float))
    up = np.array([prm.p * float(M(prm.with_(L=l), 0.0)) / prm.T0 for l in L])
    lo = np.array([prm.p * float(M_left_QD(prm.with_(L=l))) / prm.T0 for l in L])
    return up, lo
