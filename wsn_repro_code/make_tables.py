"""Regenerates Table 3 (core means), Table 5 (stress HND means + SD) and the paired Wilcoxon/Holm/CI table from outputs/outcomes_all.csv.
Test: one-sided 'greater' for FND, HND, mean_lifetime, Jain; two-sided for LND. scipy.stats.wilcoxon, zero_method='wilcox' (zero differences dropped),
normal approximation without continuity correction (method='approx', correction=False). p-values are asymptotic, not exact."""
import numpy as np, pandas as pd
from scipy import stats
O = 'outputs/'; o = pd.read_csv(O+'outcomes_all.csv'); core = o[o.phase == 'core']
ORDER = ['LEACH', 'MH-HEED', 'LEACH-GSO', 'LEACH-PSO', 'MH-HEED-GSO', 'MH-HEED-PSO', 'APSO-inspired', 'WAOA-inspired']
cols = ['FND', 'HND', 'LND', 'Jain', 'mean_lifetime', 'runtime_s']
t3 = core.groupby('protocol')[cols].mean().loc[ORDER]; t3.to_csv(O+'Table3.csv')
ci = core.groupby('protocol')[cols[:5]].agg(lambda x: stats.t.ppf(.975, len(x)-1)*x.std(ddof=1)/np.sqrt(len(x))).loc[ORDER]; ci.columns = [c+'_ci95_halfwidth' for c in ci.columns]; ci.to_csv(O+'Table3_ci95.csv')
st = o[o.phase == 'stress']; P = ['LEACH', 'MH-HEED', 'MH-HEED-PSO', 'WAOA-inspired']
g = st.groupby(['scenario', 'protocol']).HND; m, sd = g.mean().unstack()[P], g.std().unstack()[P]
m.loc[['small_low', 'large_high', 'heterogeneous', 'central_sink', 'two_cluster']].to_csv(O+'Table5.csv'); sd.loc[m.index].to_csv(O+'Table5_sd.csv')
def holm(p):
    k = len(p); r = np.argsort(p); a = np.empty(k); run = 0
    for i, j in enumerate(r): run = max(run, (k-i)*p[j]); a[j] = min(1, run)
    return a
pairs = [('LEACH-GSO', 'LEACH'), ('LEACH-PSO', 'LEACH'), ('MH-HEED-GSO', 'MH-HEED'), ('MH-HEED-PSO', 'MH-HEED'), ('MH-HEED-PSO', 'MH-HEED-GSO'), ('LEACH-PSO', 'LEACH-GSO')]
rows = []
for metric in ['FND', 'HND', 'LND', 'mean_lifetime', 'Jain']:
    for a, b in pairs:
        x = core[core.protocol == a].set_index('seed')[metric].sort_index(); y = core[core.protocol == b].set_index('seed')[metric].sort_index(); d = x-y
        alt = 'two-sided' if metric == 'LND' else 'greater'; w = stats.wilcoxon(d, alternative=alt, zero_method='wilcox', method='approx', correction=False)
        hw = stats.t.ppf(.975, 19)*d.std(ddof=1)/np.sqrt(20)
        rows.append(dict(metric=metric, optimized=a, baseline=b, alternative=alt, W=w.statistic, p_raw=w.pvalue, zero_pairs=int((d == 0).sum()), n=20, mean=d.mean(), sd=d.std(ddof=1),
                         median=d.median(), ci95_low=d.mean()-hw, ci95_high=d.mean()+hw, mean_gain_percent_of_baseline_mean=100*d.mean()/y.mean()))
S = pd.DataFrame(rows); S['p_Holm_within_metric'] = S.groupby('metric').p_raw.transform(lambda x: holm(x.values)); S['p_Holm_all_metrics'] = holm(S.p_raw.values)
S.to_csv(O+'paired_Wilcoxon_Holm_CI.csv', index=False); print(t3.round(3)); print(m.round(1))
