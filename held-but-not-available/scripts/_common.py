"""Helpers shared by the reproduction scripts."""
from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

if not os.environ.get("DISPLAY") and "ipykernel" not in sys.modules:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUTDIR = ROOT / "output"

BLUE = "#1f4e9c"
RED = "#b22222"

plt.rcParams.update({
    "font.size": 10,
    "axes.grid": True,
    "grid.color": "#e3e3e3",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
})


def outdir(path=None) -> Path:
    p = Path(path) if path is not None else OUTDIR
    p.mkdir(parents=True, exist_ok=True)
    return p


def finish(fig, name: str, path=None, show: bool = False) -> Path:
    """Save the figure to output/ (or ``path``), then show it or close it."""
    f = outdir(path) / name
    fig.savefig(f)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return f


def bisect(f, lo: float, hi: float, tol: float = 1e-10) -> float:
    """Root of a monotone function f on [lo, hi] where f(lo) and f(hi) differ in sign."""
    flo = f(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)
