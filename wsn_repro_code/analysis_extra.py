"""Analysis of the additional experiments. Run after simulator.py (all phases) and experiments.py (ctrl, ablation, repeat)."""
import numpy as np, pandas as pd
from scipy import stats
O = 'outputs/'
base = pd.read_csv(O+'outcomes_all.csv'); core = base[base.phase == 'core']
ctrl = pd.read_csv(O+'experiment_ctrl.csv'); abl = pd.read_csv(O+'experiment_ablation.csv'); rep = pd.read_csv(O+'experiment_repeat.csv')
METRICS = ['FND', 'HND', 'LND', 'mean_lifetime']

def holm(p):
    p = np.asarray(p); k = len(p); o = np.argsort(p); adj = np.empty(k); run = 0
    for i, j in enumerate(o): run = max(run, (k-i)*p[j]); adj[j] = min(1, run)
    return adj

def paired(a, b, alt):
    d = (a-b); tc = stats.t.ppf(.975, len(d)-1); hw = tc*d.std(ddof=1)/np.sqrt(len(d))
    try: w = stats.wilcoxon(d, alternative=alt, zero_method='wilcox', method='approx', correction=False); W, p = w.statistic, w.pvalue
    except ValueError: W, p = np.nan, 1.0
    return dict(mean_diff=d.mean(), ci95_low=d.mean()-hw, ci95_high=d.mean()+hw, W=W, p_raw=p, zero_pairs=int((d == 0).sum()), wins=int((d > 0).sum()), losses=int((d < 0).sum()))

def series(df, proto, metric, key='seed'):
    return df[df.protocol == proto].set_index(key)[metric].sort_index()

# ---------- E1: controls within the MH-HEED family (same pool, same routing) ----------
ctrl = ctrl.drop(columns='protocol').rename(columns={'variant': 'protocol'})  # 'variant' holds the control name; 'protocol' is the underlying engine label
allc = pd.concat([core[core.protocol.isin(['MH-HEED', 'MH-HEED-PSO', 'MH-HEED-GSO'])], ctrl[['protocol', 'seed', *METRICS, 'Jain']]])
order = ['MH-HEED', 'MH-HEED-GREEDY', 'MH-HEED-RS', 'MH-HEED-HC', 'MH-HEED-EXH', 'MH-HEED-GSO', 'MH-HEED-PSO']
summ = allc.groupby('protocol')[METRICS+['Jain']].mean().loc[order]; summ.to_csv(O+'E1_controls_means.csv')
rows = []
for m in METRICS:
    alt = 'two-sided' if m == 'LND' else 'greater'
    for a, b in [(x, 'MH-HEED') for x in order[1:]] + [('MH-HEED-PSO', x) for x in ['MH-HEED-GREEDY', 'MH-HEED-RS', 'MH-HEED-HC', 'MH-HEED-EXH', 'MH-HEED-GSO']]:
        rows.append(dict(metric=m, optimized=a, baseline=b, alternative=alt, **paired(series(allc, a, m), series(allc, b, m), alt)))
E1 = pd.DataFrame(rows); E1['p_Holm_within_metric'] = E1.groupby('metric').p_raw.transform(lambda x: holm(x)); E1['p_Holm_all'] = holm(E1.p_raw); E1.to_csv(O+'E1_controls_paired.csv', index=False)

# ---------- E2: leave-one-term-out ablation of the fitness (MH-HEED-PSO) ----------
full = core[core.protocol == 'MH-HEED-PSO'].set_index('seed'); rows = []
for v, g in abl.groupby('variant'):
    g = g.set_index('seed').sort_index()
    for m in METRICS: rows.append(dict(variant=v, metric=m, mean_variant=g[m].mean(), mean_full=full[m].mean(), **paired(g[m], full.loc[g.index, m], 'two-sided')))
E2 = pd.DataFrame(rows); E2['p_Holm_within_metric'] = E2.groupby('metric').p_raw.transform(lambda x: holm(x)); E2['p_Holm_all'] = holm(E2.p_raw); E2.to_csv(O+'E2_ablation.csv', index=False)

# ---------- E3: repeated optimizer draws per topology (10 draws, rep 0 = original) ----------
r0 = pd.concat([core[core.protocol.isin(['MH-HEED', 'MH-HEED-PSO', 'MH-HEED-GSO'])][['protocol', 'seed', *METRICS]].assign(rep=0),
                ctrl[ctrl.protocol == 'MH-HEED-RS'][['protocol', 'seed', *METRICS]].assign(rep=0)])
R = pd.concat([r0, rep.drop(columns=['protocol', 'seed']).rename(columns={'variant': 'protocol', 'base_seed': 'seed'})[['protocol', 'seed', 'rep', *METRICS]]])
assert R.groupby(['protocol', 'seed']).size().eq(10).all()
rows = []
for p, g in R.groupby('protocol'):
    for m in ['HND', 'mean_lifetime', 'FND']:
        within = g.groupby('seed')[m].std(ddof=1).pow(2).mean()**.5; between = g.groupby('seed')[m].mean().std(ddof=1)
        rows.append(dict(protocol=p, metric=m, mean=g[m].mean(), sd_within_topology_draws=within, sd_between_topology_means=between))
E3a = pd.DataFrame(rows); E3a.to_csv(O+'E3_variance_decomposition.csv', index=False)
tm = R.groupby(['protocol', 'seed'])[METRICS].mean().reset_index(); rows = []
for m in METRICS:
    alt = 'two-sided' if m == 'LND' else 'greater'
    for a, b in [('MH-HEED-PSO', 'MH-HEED'), ('MH-HEED-GSO', 'MH-HEED'), ('MH-HEED-RS', 'MH-HEED'), ('MH-HEED-PSO', 'MH-HEED-GSO'), ('MH-HEED-PSO', 'MH-HEED-RS')]:
        rows.append(dict(metric=m, optimized=a, baseline=b, alternative=alt, **paired(series(tm, a, m), series(tm, b, m), alt)))
E3b = pd.DataFrame(rows); E3b['p_Holm_within_metric'] = E3b.groupby('metric').p_raw.transform(lambda x: holm(x)); E3b['p_Holm_all'] = holm(E3b.p_raw); E3b.to_csv(O+'E3_topology_mean_paired.csv', index=False)

# ---------- E4: coverage-based service metric from node death rounds (core runs) ----------
d = pd.read_csv(O+'node_death_rounds.csv')
d = d[(d.phase == 'core')]; gx, gy = np.meshgrid(np.arange(.5, 100, 1), np.arange(.5, 100, 1)); cells = np.c_[gx.ravel(), gy.ravel()]; H = 1000; rows = []
for (sd, pr), g in d.groupby(['seed', 'protocol']):
    g = g.sort_values('node_id'); xy = g[['x_m', 'y_m']].values; dr = g.death_round.values
    dist = np.linalg.norm(cells[:, None, :]-xy[None], axis=2)
    for rs in (10, 15, 20):
        end = np.where(dist <= rs, dr[None], 0).max(axis=1)            # round until which each cell is covered
        c0 = (end > 0).mean(); t = np.arange(1, 2500)
        cov = np.array([(end >= x).mean() for x in t[::10]])
        rows.append(dict(seed=sd, protocol=pr, rs=rs, cov_initial=c0, cov_auc_H=np.minimum(end, H).mean()/H,
                         cov_life_90pct_of_initial=int(((cov >= .9*c0).sum())*10), alive_auc_H=np.minimum(dr, H).mean()/H))
E4 = pd.DataFrame(rows); E4.to_csv(O+'E4_coverage_per_run.csv', index=False)
E4s = E4.groupby(['rs', 'protocol'])[['cov_initial', 'cov_auc_H', 'cov_life_90pct_of_initial', 'alive_auc_H']].mean().reset_index(); E4s.to_csv(O+'E4_coverage_means.csv', index=False)
rows = []
for rs in (10, 15, 20):
    for m in ['cov_auc_H', 'cov_life_90pct_of_initial']:
        for a, b in [('MH-HEED-PSO', 'MH-HEED'), ('MH-HEED-GSO', 'MH-HEED'), ('MH-HEED-PSO', 'MH-HEED-GSO'), ('WAOA-inspired', 'MH-HEED-PSO'), ('MH-HEED-PSO', 'LEACH')]:
            x = E4[(E4.rs == rs)&(E4.protocol == a)].set_index('seed')[m].sort_index(); y = E4[(E4.rs == rs)&(E4.protocol == b)].set_index('seed')[m].sort_index()
            rows.append(dict(rs=rs, metric=m, optimized=a, baseline=b, **paired(x, y, 'two-sided')))
E4p = pd.DataFrame(rows); E4p['p_Holm_all'] = holm(E4p.p_raw); E4p.to_csv(O+'E4_coverage_paired.csv', index=False)

# ---------- E5: relay / receive burden from core logs ----------
rows = []
for pr in ['LEACH', 'MH-HEED', 'LEACH-GSO', 'LEACH-PSO', 'MH-HEED-GSO', 'MH-HEED-PSO', 'APSO-inspired', 'WAOA-inspired']:
    for sd in range(1001, 1021):
        z = np.load(f'{O}logs/core/base/{pr}/{sd}.npz'); span = z['span'][:, None]; rx = (z['rx_bits']*span).sum(axis=0).astype(float); role = z['role']
        ch_rounds = ((role >= 2)*span).sum(axis=0).astype(float); relay_rounds = ((role == 3)*span).sum(axis=0).astype(float)
        rows.append(dict(protocol=pr, seed=sd, rx_cv=rx.std()/rx.mean(), rx_max_share=rx.max()/rx.sum(), ch_rounds_cv=ch_rounds.std()/ch_rounds.mean(),
                         relay_rounds_total=relay_rounds.sum(), relay_nodes_used=int((relay_rounds > 0).sum())))
E5 = pd.DataFrame(rows); E5.groupby('protocol')[['rx_cv', 'rx_max_share', 'ch_rounds_cv', 'relay_rounds_total', 'relay_nodes_used']].mean().to_csv(O+'E5_relay_burden_means.csv'); E5.to_csv(O+'E5_relay_burden_per_run.csv', index=False)
print("done")
