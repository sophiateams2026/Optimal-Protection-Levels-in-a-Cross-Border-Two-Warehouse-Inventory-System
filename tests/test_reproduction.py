"""Run with:  pytest -q

The fast tests check the core formulas. The last test runs the full
reproduction and fails if any number in paper_values.json is not reproduced.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import hbna  # noqa: E402
from hbna import model, simulation  # noqa: E402


def test_table1_derivations():
    hbna.baseline().check_table1()


@pytest.mark.parametrize("R, J", [(0, 297.00), (268, 236.55), (400, 304.20)])
def test_table2_totals(R, J):
    assert model.J(hbna.baseline(), R) == pytest.approx(J, abs=0.005)


@pytest.mark.parametrize("L", [1.0, 3.0, 4.5, 8.0])
@pytest.mark.parametrize("R", [0.0, 150.0, 268.0, 400.0])
def test_closed_forms_equal_lemma2(L, R):
    prm = hbna.baseline().with_(L=L)
    a = model.expected_costs_numeric(prm, R, n=100_000)
    b = model.expected_costs_closed(prm, R)
    assert np.allclose(a, b, atol=1e-6)


@pytest.mark.parametrize("R", [0.0, 268.0, 400.0])
def test_simulation_equals_closed_forms(R):
    prm = hbna.baseline()
    a = simulation.expected_costs_sim(prm, R, n_onsets=400, dt=0.002)
    b = model.expected_costs_closed(prm, R)
    assert np.allclose(a, b, rtol=2e-3, atol=0.05)


def test_theorem1_matches_bruteforce():
    prm = hbna.baseline()
    assert model.Rstar_theorem(prm) == pytest.approx(model.Rstar_bruteforce(prm), abs=0.5)


def test_full_reproduction(tmp_path):
    import reproduce_all
    assert reproduce_all.main(["--quick", "--outdir", str(tmp_path)]) == 0


# ---------------------------------------------------------------------------
# Random parameter sets: the results must hold beyond the baseline
# ---------------------------------------------------------------------------
def _random_params(rng):
    Q = rng.uniform(200, 3000)
    d = rng.uniform(10, 500)
    return hbna.Params(
        Q=Q, d=d, QD=rng.uniform(0.05, 0.95) * Q, L=rng.uniform(0.05, 1.5) * Q / d,
        sD=rng.uniform(0, 1), sO=rng.uniform(0, 1), compliance=rng.uniform(0, 1),
        delta=rng.uniform(0, 0.01), penalty=rng.uniform(0, 5), kappa=rng.uniform(0, 20),
        lam=rng.uniform(0, 0.1) * d / Q)


@pytest.mark.parametrize("seed", range(40))
def test_random_parameters(seed):
    """Closed forms = Lemma 2; Eq. (15) = brute force; J convex; regime consistent."""
    rng = np.random.default_rng(seed)
    prm = _random_params(rng)
    for R in rng.uniform(0, prm.QD, 3):
        a = np.array(model.expected_costs_numeric(prm, R, n=50_000))
        b = np.array(model.expected_costs_closed(prm, R))
        assert np.allclose(a, b, rtol=1e-6, atol=1e-6)
    grid = np.linspace(0, prm.QD, 801)
    t = model.Rstar_theorem(prm, step=prm.QD / 4000)
    bf = model.Rstar_bruteforce(prm, grid, n=3000)
    assert abs(t - bf) <= prm.QD / 800 + 1e-9
    J = model.J_curve(prm, grid)
    assert np.diff(J, 2).min() >= -1e-6 * max(1.0, np.abs(J).max())
    reg = model.regime(prm)
    assert (reg == "DF") == (t == 0)


@pytest.mark.parametrize("seed", range(5))
def test_random_parameters_simulation(seed):
    """The time-stepping simulation agrees with the closed forms for random inputs."""
    rng = np.random.default_rng(100 + seed)
    prm = _random_params(rng)
    R = rng.uniform(0, prm.QD)
    taus = simulation.onset_grid(prm, 400)
    sim = simulation.simulate(prm, R, taus, dt=prm.T0 / 5000)
    a = np.array([sim["ES"].mean(), sim["AB"].mean(), sim["U"].mean()])
    b = np.array(model.expected_costs_on_onsets(prm, R, taus))
    w = np.array([prm.p * prm.ce, prm.p * prm.cb, prm.p * prm.kappa])
    scale = max(abs(prm.pi * R) + w @ np.abs(b), 1.0)
    assert abs(w @ (a - b)) <= 1e-3 * scale


def test_explore_cli_runs(tmp_path, monkeypatch):
    import explore
    monkeypatch.setattr(explore, "OUT", tmp_path)
    assert explore.main(["--set", "L=4", "p=0.25", "--R", "100"]) == 0
    assert explore.main(["--set", "QD=2000"]) == 2            # invalid input is rejected
    header, rows = explore.sweep(hbna.baseline(), "L", 1, 6, 4, quiet=True)
    assert len(rows) == 4
