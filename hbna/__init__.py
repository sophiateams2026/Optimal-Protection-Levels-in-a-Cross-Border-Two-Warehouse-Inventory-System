"""Replication code for "Held but Not Available: Optimal Protection Levels in a
Cross-Border Two-Warehouse Inventory System".

Modules
-------
params      Baseline parameters of Table 1 and their accounting decomposition.
model       Analytical results: Lemma 2 (state at onset), Lemma 3 (budgeted
            holding cost), Lemma 4 with the closed forms of Appendix A,
            Eq. (10) psi, Eq. (13) M(R) and Theorem 1.
simulation  An independent time-stepping simulation of the issuing policy in
            Definition 1. It produces Fig. 1 and cross-checks the formulas.
"""
from .params import Params, baseline, exporter
from . import model, simulation

__all__ = ["Params", "baseline", "exporter", "model", "simulation"]
__version__ = "1.0.0"
