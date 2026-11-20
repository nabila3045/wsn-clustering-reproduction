"""Checks (1) outputs/outcomes_all.csv equals reference/outcomes_all_reference.csv for all 416 runs (runtime/log columns ignored);
(2) outputs/node_death_rounds.csv reproduces FND/HND/LND/mean/Jain of every run; (3) no run is censored at max_rounds."""
import numpy as np, pandas as pd, sys
a = pd.read_csv('outputs/outcomes_all.csv'); b = pd.read_csv('reference/outcomes_all_reference.csv'); K = ['phase', 'scenario', 'seed', 'protocol']; C = ['nodes', 'FND', 'HND', 'LND', 'mean_lifetime', 'Jain', 'censored_nodes']
m = a.merge(b, on=K, suffixes=('_n', '_r')); ok = len(a) == len(b) == len(m) == 416
for c in C: ok &= bool(np.allclose(m[c+'_n'], m[c+'_r'], rtol=0, atol=1e-9))
print('outcomes equal to reference:', ok)
d = pd.read_csv('outputs/node_death_rounds.csv'); bad = 0
for k, g in d.groupby(K):
    r = a[(a.phase == k[0]) & (a.scenario == k[1]) & (a.seed == k[2]) & (a.protocol == k[3])].iloc[0]; t = np.sort(g.death_round.values); n = len(t); h = (n+1)//2
    bad += not (t[0] == r.FND and t[h-1] == r.HND and t[-1] == r.LND and abs(t.mean()-r.mean_lifetime) < 1e-6 and abs(t.sum()**2/(n*(t*t).sum())-r.Jain) < 1e-9 and (g.censored == 0).all())
print('death-round file mismatches:', bad, '| censored runs:', int((a.censored_nodes > 0).sum()))
sys.exit(0 if ok and bad == 0 else 1)
