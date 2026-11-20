# WSN clustering simulator and evaluation — reproduction package

Code, configuration, seed schedule, per-seed outcomes, per-node death rounds, statistics and plotting scripts for the manuscript
"Protocol Driven PSO and GSO Optimization of LEACH and Multi Hop HEED ..." (IJIES).

## Status of this simulator (read first)
`simulator.py` is an **independent reconstruction** with fully specified rules (`MODEL_SPEC.md`). It is not the historical source code, which is not available.
Every detail not specified earlier was fixed explicitly in this reconstruction. Results produced here therefore **supersede** any numbers obtained with the earlier code.
All 416 runs reproduce exactly (`verify_reproduction.py`). The model is simulation only: first-order radio model, no radio range cap, no collisions, no retransmissions,
centralised state collection/dissemination cost not modelled, round-end energy debit (no partial-packet delivery). It supports no claim about hardware or packet delivery.
`APSO-inspired` and `WAOA-inspired` are local approximations of the cited methods, not replications.

## Run
```
pip install -r requirements.txt      # python 3.13
./run_all.sh                         # ~6 min on 2 cores; deterministic (seeds in config.json, stream map in seed_stream_mapping.csv)
```
Single steps: `python3 simulator.py --phase core|stress|development|diagnostic|all --workers N`.

## Layout
| File | Content |
|---|---|
| `simulator.py`, `config.json` | simulator and all parameters (Table 2 values, fitness weights, PSO/GSO coefficients, stress scenarios, seeds) |
| `MODEL_SPEC.md` | executable specification, hand-checked example, log format |
| `seed_stream_mapping.csv` | random-stream entropy per phase/scenario/seed/protocol (topology and energy streams are shared by all protocols of a seed) |
| `outputs/outcomes_all.csv` | one row per run (416): FND, HND, LND, mean lifetime, Jain, runtime |
| `outputs/node_death_rounds.csv` | one row per node and run (21,800): coordinates, death round, CH appointments |
| `outputs/logs/**.npz` | per-run event-block logs (energy, parent, role, TX/RX/aggregation bits) — archive these in the DOI release |
| `make_tables.py` | Table 3, Table 5 (+SD), paired Wilcoxon/Holm/CI (`outputs/paired_Wilcoxon_Holm_CI.csv`) |
| `experiments.py`, `analysis_extra.py` | added controls: `MH-HEED-GREEDY`, `-HC`, `-RS`, `-EXH`; leave-one-term-out ablation; 10 optimizer draws per topology; coverage and relay metrics (`outputs/E1…E5*.csv`) |
| `analysis_sens.py` | stress-scenario paired contrasts (`outputs/stress_paired.csv`) and one-at-a-time sensitivity (sink, payload, interval; `experiments.py sens` → `outputs/sens_*.csv`) |
| `benchmark_runtime.py` | runtime of all MH-HEED variants under identical conditions (`outputs/E6_runtime.csv`) |
| `verify_reproduction.py`, `tests/test_hand_example.py` | consistency checks |
| `make_figures.py` | figures (`outputs/figures`) |

## Statistics
Paired by topology seed. Wilcoxon signed-rank (`scipy.stats.wilcoxon`, zero differences dropped, normal approximation without continuity correction — p-values are asymptotic, not exact).
One-sided ("greater") for FND, HND, mean lifetime, Jain; two-sided for LND; for ablations and coverage metrics two-sided. Holm adjustment within metric and across all tests of a table.
Confidence intervals: two-sided Student-t, 95%, on paired differences. No run reached the 5,000-round limit (largest death round 2,493), so nothing is censored.
Minor note: one p-value (mean lifetime, LEACH-PSO vs LEACH-GSO, 0.05218 vs 0.05216 in an earlier table) depends on the scipy version.

## Coverage metric (assumption)
Coverage uses an explicit sensing radius r_s (10, 15, 20 m; main text value 15 m) on a 1 m grid of the 100x100 m field; a cell is covered until its last covering node dies.
Because the model has no radio range cap and re-elects CHs immediately after a CH death, every alive node always has a route to the sink; the "sink-reachable alive fraction" equals the alive fraction.
