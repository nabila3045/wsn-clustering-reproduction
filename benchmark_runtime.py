"""Software runtime of the whole simulation per protocol under identical conditions (single process, no log writing), core seeds 1001-1020.
Wall-clock simulation time on the host, including objective evaluations; it is NOT on-node computation time or energy."""
import time, numpy as np, pandas as pd
import experiments as x, simulator as s
P = ['MH-HEED', 'MH-HEED-GREEDY', 'MH-HEED-RS', 'MH-HEED-HC', 'MH-HEED-EXH', 'MH-HEED-GSO', 'MH-HEED-PSO']; rows = []
for rep in range(3):
    for p in P:
        for sd in s.CFG['core_seeds']:
            t = time.perf_counter(); s.simulate(('bench', 'base', sd, p)); rows.append(dict(protocol=p, rep=rep, seed=sd, seconds=time.perf_counter()-t))
d = pd.DataFrame(rows); g = d.groupby('protocol').seconds.agg(['mean', 'std']).loc[P]; g['ratio_vs_MH-HEED'] = g['mean']/g.loc['MH-HEED', 'mean']; g.to_csv('outputs/E6_runtime.csv'); print(g.round(4))
