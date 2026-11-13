"""Table 1: baseline parameters and their accounting derivation.

Checks that every derived entry of Table 1 (T0, QO, Delta s, c_e, c_b, p)
follows from its components, and writes output/table1.csv.

Run:  python scripts/table1_parameters.py
"""
from __future__ import annotations

import csv

from _common import outdir
import hbna


def run(path=None) -> dict:
    prm = hbna.baseline()
    checks = prm.check_table1()  # raises AssertionError on any mismatch
    with open(outdir(path) / "table1.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Parameter", "Value", "Derivation"])
        w.writerows(prm.table1_rows())
    res = {"table1." + k.split(" =")[0].replace(" ", ""): v[0] for k, v in checks.items()}
    res.update({"table1.pi": prm.pi, "table1.X": prm.X})
    return res


if __name__ == "__main__":
    prm = hbna.baseline()
    for row in prm.table1_rows():
        print(f"{row[0]:14s} {row[1]:22s} {row[2]}")
    out = run()
    print("\nAll derived entries of Table 1 check out:")
    for k, v in out.items():
        print(f"  {k:24s} {v:.4f}")
