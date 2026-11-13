# Held but Not Available: replication code

Code for the manuscript *Held but Not Available: Optimal Protection Levels in a Cross-Border Two-Warehouse Inventory System*.

The code does two things. It reproduces every table, figure and number in Section 6. It also lets you change any input in Table 1 and see the result, with each answer checked by independent computations.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/sophiateams2026/Optimal-Protection-Levels-in-a-Cross-Border-Two-Warehouse-Inventory-System/blob/main/notebooks/reproduce_and_explore.ipynb)
[![reproduce](https://github.com/sophiateams2026/Optimal-Protection-Levels-in-a-Cross-Border-Two-Warehouse-Inventory-System/actions/workflows/reproduce.yml/badge.svg)](https://github.com/sophiateams2026/Optimal-Protection-Levels-in-a-Cross-Border-Two-Warehouse-Inventory-System/actions/workflows/reproduce.yml)

## Quick start

Python 3.9 or later. Only NumPy and Matplotlib are needed.

```bash
git clone https://github.com/sophiateams2026/Optimal-Protection-Levels-in-a-Cross-Border-Two-Warehouse-Inventory-System.git
cd Optimal-Protection-Levels-in-a-Cross-Border-Two-Warehouse-Inventory-System
pip install -r requirements.txt
python reproduce_all.py          # about 20 seconds
```

The last line prints `71 of 71 checks passed.` Results go to `output/`. The file `output/verification_report.md` lists every number reported in the manuscript next to the value computed here.

No installation? Click **Open in Colab** above and choose *Runtime → Run all*.

## Try your own numbers

`explore.py` re-solves the model for any parameter values.

```bash
python explore.py --list                       # parameter names and baseline values
python explore.py --set L=4 lam=0.03           # longer and more frequent restrictions
python explore.py --set p=0.3 kappa=10         # set p directly instead of lam * T0
python explore.py --set dv=5                   # exporter case of Section 6.4
python explore.py --set ce=0.10 cb=2.5         # set cost rates directly
python explore.py --R 150 300                  # add reserves of your choice to the cost table
python explore.py --sweep L 0.5 8 31           # R* and costs over a grid of L
python explore.py --config my_params.json      # copy and edit params_baseline.json
```

For each parameter set the script prints:

* the regime (DF, partial reserve or OF) and R\* from Theorem 1, Eq. (15);
* R\* from a brute-force search over R that uses Lemma 2 only;
* a cost breakdown in the layout of Table 2;
* a cross-check of the expected costs by three routes (below), with PASS or CHECK FAILED.

It also saves the curve J(R) − J(0) to `output/explore/`. Values outside the model are rejected with a message (for example QD ≥ Q, or p > 1). Values that stretch an assumption trigger a note (for example p > 0.3, or L longer than the cycle).

The Colab notebook has the same tool in an editable cell, plus sliders.

### Parameter names

| Name | Baseline | Meaning |
|---|---:|---|
| `Q`, `d`, `QD` | 1000, 100, 400 | replenishment quantity, weekly demand, stock placed in the demand-side warehouse |
| `v`, `P`, `r` | 50, 80, 0.002 | landed unit cost, price, cost of capital per week |
| `dv` | 0 | carrying-amount gap v_D − v_O (0 for an importer, 5 for the exporter) |
| `sD`, `sO` | 0.20, 0.08 | storage rates per unit-week |
| `compliance`, `delta` | 0.10, 0.001 | c_e = compliance + delta × v |
| `penalty` | 1.50 | c_b = r v + r P + penalty |
| `kappa` | 5 | per-unit backlog cost |
| `lam`, `L` | 0.02, 3 | restriction rate per week (p = lam × T0), restriction duration in weeks |
| `p`, `ce`, `cb`, `dh` | — | optional shortcuts that set these values directly |

## How the results are verified

Expected costs are computed by three independent routes.

| Route | What it does | Code |
|---|---|---|
| A | Lemma 2 applied to each onset τ, averaged over 200,000 onsets in [0, T0) | `model.expected_costs_numeric` |
| B | Closed forms of Lemma 4 (Appendix A) | `model.expected_costs_closed` |
| C | Time-stepping simulation of the issuing rule in Definition 1. It uses no formula from the paper. | `simulation.simulate` |

Route C moves both warehouses forward in small time steps. It applies rules (i)–(iii) of Definition 1 and Assumption 2 literally: stranded stock is frozen, backlog is served first when access returns, and the next replenishment clears any remaining backlog. Fig. 1 is drawn from this simulation.

On top of this, the code checks the analytical results directly:

* R\* from Eq. (15) against brute-force minimisation (Theorem 1);
* convexity of J(R) and the slope formula Eq. (14);
* the budgeted holding cost Eq. (7) against the simulated holding cost (Lemma 3);
* the regime map of Fig. 4 against brute force on a 40 × 40 grid;
* monotonicity of R\* (Proposition 2), and the shape of V(p) (Proposition 3);
* the critical fractile of Corollary 3 with a random restriction duration.

`pytest -q` runs the same checks on 40 random parameter sets, including restrictions longer than the cycle.

## What reproduces what

| Manuscript | Script | Output |
|---|---|---|
| Table 1 | `scripts/table1_parameters.py` | `table1.csv` |
| Fig. 1, Section 6.2 | `scripts/fig1_accessible_inventory.py` | `fig1.png` |
| Table 2, Fig. 2, Section 6.3 | `scripts/table2_fig2_expected_cost.py` | `table2.csv`, `fig2.png` |
| Fig. 3, Section 6.4 | `scripts/fig3_sensitivity.py` | `fig3a.png`, `fig3b.png` |
| Fig. 4 and exporter scenario, Section 6.4 | `scripts/fig4_regimes_exporter.py` | `fig4.png` |
| Supplementary: exporter with a smaller c_e (not in the manuscript) | `scripts/supp_exporter_ce.py` | `supp_exporter_ce.csv`, `.png` |

Each script runs on its own, for example `python scripts/fig3_sensitivity.py`.

## Code and manuscript

| Manuscript | Function in `hbna/model.py` |
|---|---|
| Lemma 2, Eqs. (5)–(6) | `onset_state`, `pathwise_outcome` |
| Lemma 3, Eqs. (7)–(8) | `budgeted_holding_cost`, `Params.pi` |
| Eq. (9) | `J`, `J_components` |
| Eq. (10) | `psi` |
| Lemma 4, Eq. (11), Appendix A | `expected_E0`, `expected_costs_closed` |
| Eq. (13) | `M`, `M_left_QD` |
| Theorem 1, Eqs. (14)–(15) | `J_right_derivative`, `Rstar_theorem`, `regime` |
| Table 1 | `hbna/params.py` |

## Notes

* Table 2 reports R\* = 268, the best reserve on the integer grid. Eq. (15) gives 267.90 on a continuous scale. The difference in J is below 0.001.
* L is deterministic, as in the baseline of Section 6. Corollary 3 is checked separately with a random L.
* J(R) is defined net of terms that do not depend on R. It can be zero or negative (for example when p = 0 or Δh < 0). Percentage savings are then not reported.
* Route C has an error of order dt from rounding event times to the time grid. In the baseline the onsets fall on the grid and the error is below 10⁻¹⁰.

## Repository layout

```
hbna/                 model code (params.py, model.py, simulation.py)
scripts/              one script per table or figure
reproduce_all.py      runs every script and writes the verification report
explore.py            try your own parameter values
paper_values.json     numbers reported in the manuscript, with tolerances
params_baseline.json  Table 1 as an editable template
notebooks/            Colab notebook
tests/                pytest suite
output/               results produced by reproduce_all.py
```

## Citation and licence

See `CITATION.cff`. Code released under the MIT licence.
