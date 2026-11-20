"""Stress-scenario paired contrasts (Table 5 companion) and one-at-a-time sensitivity (outputs/experiment_sens.csv)."""
import numpy as np, pandas as pd
from scipy import stats
O = 'outputs/'
def holm(p):
    p = np.asarray(p); k = len(p); o = np.argsort(p); a = np.empty(k); run = 0
    for i, j in enumerate(o): run = max(run, (k-i)*p[j]); a[j] = min(1, run)
    return a
def paired(x, y):
    d = x-y; hw = stats.t.ppf(.975, len(d)-1)*d.std(ddof=1)/np.sqrt(len(d))
    try: p = stats.wilcoxon(d, alternative='two-sided', zero_method='wilcox', method='approx', correction=False).pvalue
    except ValueError: p = 1.0
    return dict(mean_diff=d.mean(), ci95_low=d.mean()-hw, ci95_high=d.mean()+hw, p_raw=p, wins=int((d > 0).sum()), losses=int((d < 0).sum()), ties=int((d == 0).sum()))
o = pd.read_csv(O+'outcomes_all.csv'); st = o[o.phase == 'stress']; rows = []
for sc, g in st.groupby('scenario', sort=False):
    w = g.pivot(index='seed', columns='protocol', values='HND')
    for a, b in [('MH-HEED-PSO', 'MH-HEED'), ('MH-HEED-PSO', 'LEACH'), ('WAOA-inspired', 'MH-HEED-PSO')]: rows.append(dict(scenario=sc, optimized=a, baseline=b, **paired(w[a], w[b])))
S = pd.DataFrame(rows); S['p_Holm_all'] = holm(S.p_raw); S.to_csv(O+'stress_paired.csv', index=False)
s = pd.read_csv(O+'experiment_sens.csv').drop(columns='protocol').rename(columns={'variant': 'setting'}); s['protocol'] = s.merge(pd.read_csv(O+'experiment_sens.csv')[['variant']], left_index=True, right_index=True)['variant'] if False else pd.read_csv(O+'experiment_sens.csv')['protocol']
rows = []; means = []
for sc, g in s.groupby('setting', sort=False):
    for m in ['HND', 'mean_lifetime']:
        w = g.pivot(index='seed', columns='protocol', values=m); means.append(dict(setting=sc, metric=m, **{p: w[p].mean() for p in w.columns}))
        for a, b in [('MH-HEED-PSO', 'MH-HEED'), ('MH-HEED-RS', 'MH-HEED'), ('MH-HEED-PSO', 'MH-HEED-RS'), ('MH-HEED-PSO', 'LEACH')]: rows.append(dict(setting=sc, metric=m, optimized=a, baseline=b, **paired(w[a], w[b])))
P = pd.DataFrame(rows); P['p_Holm_all'] = holm(P.p_raw); P.to_csv(O+'sens_paired.csv', index=False); pd.DataFrame(means).to_csv(O+'sens_means.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.float_format', lambda x: f'{x:.4g}')
print(pd.DataFrame(means)); print(P[(P.metric == 'HND')][['setting', 'optimized', 'baseline', 'mean_diff', 'ci95_low', 'ci95_high', 'wins', 'losses', 'ties', 'p_raw', 'p_Holm_all']]); print(S)
