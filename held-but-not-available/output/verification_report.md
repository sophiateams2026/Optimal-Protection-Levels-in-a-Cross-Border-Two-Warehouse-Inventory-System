# Verification report

Generated 2026-11-13 13:26 UTC with Python 3.12.3, NumPy 2.5.4, Matplotlib 3.11.2, hbna 1.0.0.

**71 of 71 checks passed.**

## A. Numbers reported in the manuscript

| Where | Statement | Reported | Computed | Tolerance | Result |
|---|---|---:|---:|---:|:---:|
| Table 1 | T0 = Q/d = 10 weeks | 10 | 10 | 1e-09 | PASS |
| Table 1 | QO = Q - QD = 600 | 600 | 600 | 1e-09 | PASS |
| Table 1 | Delta s = 0.20 - 0.08 = 0.12 | 0.12 | 0.12 | 1e-09 | PASS |
| Table 1 | c_e = 0.10 + 0.001 x 50 | 0.15 | 0.15 | 1e-09 | PASS |
| Table 1 | c_b = rv 0.10 + rP 0.16 + penalty 1.50 | 1.76 | 1.76 | 1e-09 | PASS |
| Table 1 | p = lambda T0 = 0.02 x 10 | 0.2 | 0.2 | 1e-09 | PASS |
| Table 1 | r = 0.2%/wk is about 10.4% a year | 0.104 | 0.104 | 1e-09 | PASS |
| Section 6.2, Fig. 1 | Accessible inventory drops to 200 under DF at the onset | 200 | 200 | 0.5 | PASS |
| Section 6.2, Fig. 1 | ... 300 under P(300) | 300 | 300 | 0.5 | PASS |
| Section 6.2, Fig. 1 | ... 400 under OF | 400 | 400 | 0.5 | PASS |
| Section 6.2, Fig. 1 | DF runs out at week 4 | 4 | 4 | 0.01 | PASS |
| Section 6.2, Fig. 1 | DF backlogs 100 units before access returns at week 5 | 100 | 100 | 0.5 | PASS |
| Section 6.2, Fig. 1 | P(300) just covers the restriction | 0 | 0 | 0.5 | PASS |
| Section 6.2, Fig. 1 | OF still holds 100 demand-side units at restoration | 100 | 100.1 | 0.5 | PASS |
| Section 6.2, Fig. 1 | Minimum physical stock during the restriction over the three policies (DF 600; P(300) and OF 500) | 500 | 500.1 | 0.5 | PASS |
| Table 2 | R* = 268 (integer grid) | 268 | 268 | 0 | PASS |
| Table 2 | DF resilience premium | 0 | 0 | 0.005 | PASS |
| Table 2 | DF stranded-inventory cost | 37.8 | 37.8 | 0.005 | PASS |
| Table 2 | DF backlog-time cost | 79.2 | 79.2 | 0.005 | PASS |
| Table 2 | DF per-unit backlog cost | 180 | 180 | 0.005 | PASS |
| Table 2 | DF expected relevant cost | 297 | 297 | 0.005 | PASS |
| Table 2 | DF memo E[ES] | 1260 | 1260 | 0.05 | PASS |
| Table 2 | DF memo E[U] | 180 | 180 | 0.005 | PASS |
| Table 2 | P(R*) resilience premium | 192.96 | 192.96 | 0.005 | PASS |
| Table 2 | P(R*) stranded-inventory cost | 23.33 | 23.328 | 0.005 | PASS |
| Table 2 | P(R*) backlog-time cost | 1.06 | 1.0621 | 0.005 | PASS |
| Table 2 | P(R*) per-unit backlog cost | 19.2 | 19.2 | 0.005 | PASS |
| Table 2 | P(R*) expected relevant cost | 236.55 | 236.5501 | 0.005 | PASS |
| Table 2 | P(R*) memo E[ES] | 777.6 | 777.6 | 0.05 | PASS |
| Table 2 | P(R*) memo E[U] | 19.2 | 19.2 | 0.005 | PASS |
| Table 2 | OF resilience premium | 288 | 288 | 0.005 | PASS |
| Table 2 | OF stranded-inventory cost | 16.2 | 16.2 | 0.005 | PASS |
| Table 2 | OF backlog-time cost | 0 | 0 | 0.005 | PASS |
| Table 2 | OF per-unit backlog cost | 0 | 0 | 0.005 | PASS |
| Table 2 | OF expected relevant cost | 304.2 | 304.2 | 0.005 | PASS |
| Table 2 | OF memo E[ES] | 540 | 540 | 0.05 | PASS |
| Table 2 | OF memo E[U] | 0 | 0 | 0.005 | PASS |
| Section 6.3 | DF expected penalty and backlog costs are 259.20 | 259.2 | 259.2 | 0.005 | PASS |
| Section 6.3; Abstract | J(R*) is 20.4% below DF | 20.4 | 20.3535 | 0.05 | PASS |
| Section 6.3; Abstract | J(R*) is 22.2% below OF | 22.2 | 22.2386 | 0.05 | PASS |
| Section 6.3 | Comparing only OF and DF would favour DF by 7.20 | 7.2 | 7.2 | 0.005 | PASS |
| Section 6.3, Fig. 2 | At p = 0.1 DF is optimal | 0 | 0 | 0 | PASS |
| Section 6.3, Fig. 2 | At p = 0.2 the curve falls to an interior minimum | 268 | 268 | 0 | PASS |
| Section 6.3, Fig. 2 | At p = 0.3 the minimum sits at R = 300 = dL | 300 | 300 | 0 | PASS |
| Section 6.4, Fig. 3(a) | At p = 0.2 no reserve is held for restrictions shorter than about 0.55 weeks | 0.55 | 0.5465 | 0.005 | PASS |
| Section 6.4, Fig. 3(a) | R* reaches Q_D at L of about 4.2 weeks | 4.2 | 4.2129 | 0.05 | PASS |
| Section 6.4, Fig. 3(b) | For L = 3, DF is optimal up to p of about 0.128 | 0.128 | 0.1275 | 0.0005 | PASS |
| Section 6.4, Fig. 3(b) | R* = dL = 300 from p of about 0.22 onward | 0.22 | 0.2202 | 0.005 | PASS |
| Section 6.4; Corollary 2 | Full OF issuing is never optimal at L = 3 | 300 | 300 | 0 | PASS |
| Section 6.4; Corollary 2 | p c_e L = 0.45 p never exceeds Delta s T0 = 1.2 | 0.45 | 0.45 | 1e-09 | PASS |
| Section 6.4, Fig. 4 | At Delta s = 0.12 the partial region starts at L of about 0.55 | 0.55 | 0.5465 | 0.005 | PASS |
| Section 6.4, Fig. 4 | ... and ends at L of about 4.2 | 4.2 | 4.2129 | 0.05 | PASS |
| Section 6.4 | Delta h = 0.12 + 0.002 x 5 = 0.13 | 0.13 | 0.13 | 1e-09 | PASS |
| Section 6.4 | Exporter stays in the partial-protection region (1 = partial) | 1 | 1 | 0 | PASS |
| Section 6.4 | The optimal reserve falls from 268 to 237 | 237 | 237 | 0 | PASS |
| Section 6.4 | Expected relevant cost is 251.71 | 251.71 | 251.7066 | 0.005 | PASS |
| Section 6.4 | ... 15.3% below DF | 15.3 | 15.2503 | 0.05 | PASS |
| Section 6.4 | ... 23.3% below OF | 23.3 | 23.307 | 0.05 | PASS |

## B. Internal consistency checks of the analytical results

| Check | Computed | Criterion | Result |
|---|---:|---|:---:|
| Lemma 4 closed forms equal Lemma 2 integrated over the onset (L = 1, 3, 4.5, 8; 21 reserves each) | 3.54e-08 | <= 1e-06 | PASS |
| Table 2: time-stepping simulation of Definition 1 agrees with the closed forms (error of order dt) | 8.8e-05 | <= 0.01 | PASS |
| J(R) is convex on [0, Q_D] for p = 0.1, 0.2, 0.3 (Theorem 1) | -1.14e-13 | >= -1e-09 | PASS |
| Finite-difference slope of J equals Eq. (14) | 1.57e-07 | <= 0.0001 | PASS |
| R* from Eq. (15) equals brute-force minimisation of Lemma 2 costs (grid step 0.5) | 0.1 | <= 0.5 | PASS |
| Simulated holding cost of an undisrupted cycle equals Eq. (7) | 1.85e-10 | <= 1e-06 | PASS |
| Budgeted holding cost rises by pi = Delta s Q_O / d per unit of reserve (Eq. 8) | 0 | |x| <= 1e-09 | PASS |
| Regime map from Theorem 1 matches brute force on a 40 x 40 grid | 0 | <= 0 | PASS |
| R* is monotone in p, L, c_e, c_b, kappa and Delta h as in Proposition 2 | 0 | <= 0 | PASS |
| V(p) = 0 for p <= Delta s T0 / M(0) (Proposition 3) | 0 | |x| <= 1e-09 | PASS |
| V(p) is nondecreasing (Proposition 3) | 0 | >= -1e-09 | PASS |
| V(p) is convex (Proposition 3) | -2.84e-13 | >= -1e-09 | PASS |
| Corollary 3 critical fractile matches brute force when L is random (uniform on 1..5 weeks, c_b = 0) | 200 | = 200 | PASS |
