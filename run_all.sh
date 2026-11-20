#!/usr/bin/env bash
# Full regeneration (about 6 minutes on 2 cores). Run from the repository root.
set -e
python3 simulator.py --workers 2 --phase all        # 416 runs -> outputs/outcomes_all.csv, outputs/logs/**.npz
python3 export_node_death_rounds.py                  # outputs/node_death_rounds.csv
python3 verify_reproduction.py                       # compares with reference/outcomes_all_reference.csv
python3 tests/test_hand_example.py                   # hand-checkable multi-hop example (1.304 mJ)
python3 make_tables.py                               # Table 3, Table 5, paired Wilcoxon/Holm/CI
python3 experiments.py ctrl && python3 experiments.py ablation && python3 experiments.py repeat && python3 experiments.py sens && python3 analysis_sens.py
python3 analysis_extra.py                            # E1-E5 tables
python3 benchmark_runtime.py                         # E6 runtime under identical conditions
python3 make_figures.py
