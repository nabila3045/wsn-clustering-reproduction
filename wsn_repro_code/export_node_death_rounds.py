"""Writes outputs/node_death_rounds.csv (one row per node and run) from outputs/outcomes_all.csv and the per-run logs in outputs/logs."""
import numpy as np, pandas as pd
o = pd.read_csv('outputs/outcomes_all.csv'); rows = []
for r in o.itertuples():
    z = np.load(r.log); n = len(z['death']); d = z['death']
    rows.append(pd.DataFrame(dict(phase=r.phase, scenario=r.scenario, seed=r.seed, protocol=r.protocol, node_id=np.arange(n), x_m=z['xy'][:, 0], y_m=z['xy'][:, 1],
                                  initial_j=z['initial'], death_round=d, censored=(d == 0).astype(int), CH_appointments=z['ch_duty'])))
pd.concat(rows).to_csv('outputs/node_death_rounds.csv', index=False); print('written', sum(map(len, rows)), 'rows')
