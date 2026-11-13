"""Baseline parameters (Table 1 of the manuscript).

The time unit is one week. Each cost rate is stored as the accounting
components it is built from, so the derivation column of Table 1 can be
checked with :meth:`Params.check_table1`.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional


@dataclass(frozen=True)
class Params:
    # Quantities
    Q: float = 1000.0         # replenishment quantity (units)
    d: float = 100.0          # demand rate (units/week)
    QD: float = 400.0         # allocation to the demand-side warehouse (= capacity W_D)

    # Carrying amounts and price
    v: float = 50.0           # landed unit cost under IAS 2; v_O in the export case
    dv: float = 0.0           # carrying-amount gap v_D - v_O (0 for an importer)
    P: float = 80.0           # selling price
    r: float = 0.002          # cost of capital per week (about 10.4% a year)

    # Storage rates (per unit-week)
    sD: float = 0.20          # demand-side warehouse
    sO: float = 0.08          # origin-side warehouse

    # Stranded-inventory cost c_e = compliance + delta * v
    compliance: float = 0.10  # compliance, war-risk insurance, bonded fees
    delta: float = 0.001      # hazard rate of a write-down to NRV (per week)

    # Time-based backlog cost c_b = r v + r P + penalty
    penalty: float = 1.50     # contractual late-delivery penalty (per unit-week)

    kappa: float = 5.0        # per-unit backlog cost (per unit backordered)

    # Disruption process
    lam: float = 0.02         # restriction arrival rate (per week)
    L: float = 3.0            # restriction duration in weeks (deterministic baseline)

    # Optional direct overrides (None = use the accounting decomposition)
    ce_override: Optional[float] = None
    cb_override: Optional[float] = None
    p_override: Optional[float] = None
    dh_override: Optional[float] = None

    # ---- derived quantities -------------------------------------------------
    @property
    def T0(self) -> float:
        """Cycle length Q/d."""
        return self.Q / self.d

    @property
    def QO(self) -> float:
        """Overflow stored on the origin side, Q - Q_D."""
        return self.Q - self.QD

    @property
    def ds(self) -> float:
        """Storage-rate differential s_D - s_O."""
        return self.sD - self.sO

    @property
    def dh(self) -> float:
        """Holding-cost differential, Delta h = Delta s + r * Delta v (Section 3.4).

        Equals Delta s for an importer (Delta v = 0). All results of Section 5
        hold with Delta s replaced by Delta h.
        """
        if self.dh_override is not None:
            return self.dh_override
        return self.ds + self.r * self.dv

    @property
    def ce(self) -> float:
        """Stranded-inventory cost per unit-week."""
        if self.ce_override is not None:
            return self.ce_override
        return self.compliance + self.delta * self.v

    @property
    def cb(self) -> float:
        """Time-based backlog cost per unit-week."""
        if self.cb_override is not None:
            return self.cb_override
        return self.r * self.v + self.r * self.P + self.penalty

    @property
    def p(self) -> float:
        """Probability that a restriction begins in a cycle, p = lambda * T0."""
        if self.p_override is not None:
            return self.p_override
        return self.lam * self.T0

    @property
    def X(self) -> float:
        """Demand during the restriction, X = d * L."""
        return self.d * self.L

    @property
    def pi(self) -> float:
        """Resilience premium per unit of reserve, Eq. (8): Delta h * Q_O / d."""
        return self.dh * self.QO / self.d

    def with_(self, **changes) -> "Params":
        """Copy with some values changed, e.g. ``prm.with_(L=4, p=0.3)``.

        ``p``, ``ce``, ``cb`` and ``dh`` set the value directly instead of
        through their accounting components.
        """
        for key in ("p", "ce", "cb", "dh"):
            if key in changes:
                changes[f"{key}_override"] = changes.pop(key)
        return replace(self, **changes)

    # ---- Table 1 ------------------------------------------------------------
    def table1_rows(self):
        """Rows of Table 1 as (parameter, value, derivation)."""
        return [
            ("Q, d, T0", f"{self.Q:g}, {self.d:g}/wk, {self.T0:g} wk",
             "Replenishment quantity, demand rate, cycle length"),
            ("QD = WD, QO", f"{self.QD:g}, {self.QO:g}",
             "Demand-side capacity binds; overflow stored on the origin side"),
            ("v, P", f"{self.v:g}, {self.P:g}", "Landed unit cost (IAS 2) and selling price"),
            ("r", f"{100 * self.r:g}%/wk",
             f"Cost of capital, about {100 * self.r * 52:.1f}% a year"),
            ("sD, sO (ds)", f"{self.sD:.2f}, {self.sO:.2f} ({self.ds:.2f})",
             "Storage rates from warehouse cost pools"),
            ("ce", f"{self.ce:.2f}",
             f"Compliance, insurance, bonded fees {self.compliance:.2f}; "
             f"expected write-down delta*v = {self.delta:g} x {self.v:g}"),
            ("cb", f"{self.cb:.2f}",
             f"Capital rv = {self.r * self.v:.2f}; deferred revenue rP = {self.r * self.P:.2f}; "
             f"contractual penalty {self.penalty:.2f}"),
            ("kappa", f"{self.kappa:g}", "Expediting, handling and goodwill per unit backordered"),
            ("p", f"{self.p:.2f}",
             f"lambda = {self.lam:g} per week (one restriction every {1 / self.lam:g} weeks)"),
            ("L", f"{self.L:g} wk", "Restriction duration (deterministic in the baseline)"),
        ]

    def check_table1(self) -> dict:
        """Verify every derived entry of Table 1; raise AssertionError on a mismatch."""
        checks = {
            "T0 = Q/d = 10": (self.T0, 10.0),
            "QO = Q - QD = 600": (self.QO, 600.0),
            "ds = sD - sO = 0.12": (self.ds, 0.12),
            "ce = 0.10 + 0.001 x 50 = 0.15": (self.ce, 0.15),
            "cb = rv + rP + 1.50 = 1.76": (self.cb, 1.76),
            "p = lambda x T0 = 0.20": (self.p, 0.20),
            "52 r = 10.4% a year": (52 * self.r, 0.104),
        }
        for name, (got, want) in checks.items():
            if abs(got - want) > 1e-9:
                raise AssertionError(f"Table 1 check failed: {name}: got {got}")
        return checks


def baseline() -> Params:
    """Importer baseline of Table 1."""
    return Params()


def exporter() -> Params:
    """Exporter scenario of Section 6.4: v_O = 50, v_D = 55, so Delta h = 0.13."""
    return Params(dv=5.0)
