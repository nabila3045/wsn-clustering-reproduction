"""Figures from outputs/ (matplotlib). Fig1: HND distribution (core, 20 paired seeds). Fig2: optimizer convergence, diagnostic seed 777
(best-so-far fitness divided by the best of the initial population, averaged over all elections of the run; normalised per protocol, so curves are not comparable in level across protocols).
Fig3: HND of the controls within the MH-HEED family (same candidate pool and routing), mean with 95% CI."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
O = 'outputs/'; import os; os.makedirs(O+'figures', exist_ok=True)
o = pd.read_csv(O+'outcomes_all.csv'); core = o[o.phase == 'core']
ORDER = ['LEACH', 'MH-HEED', 'LEACH-GSO', 'LEACH-PSO', 'MH-HEED-GSO', 'MH-HEED-PSO', 'APSO-inspired', 'WAOA-inspired']
plt.rcParams.update({'font.family': 'serif', 'font.size': 9})
fig, ax = plt.subplots(figsize=(7, 3.4)); data = [core[core.protocol == p].HND.values for p in ORDER]
ax.boxplot(data, showmeans=True, meanprops=dict(marker='^', markerfacecolor='green', markeredgecolor='green'), medianprops=dict(color='tab:orange'))
ax.set_xticks(range(1, 9)); ax.set_xticklabels(ORDER, rotation=30, ha='right'); ax.set_ylabel('HND (rounds)'); ax.grid(alpha=.25, axis='y'); fig.tight_layout(); fig.savefig(O+'figures/fig1_hnd_boxplot.png', dpi=300)
fig, ax = plt.subplots(figsize=(7, 3.4))
for p in ['LEACH-GSO', 'LEACH-PSO', 'MH-HEED-GSO', 'MH-HEED-PSO', 'APSO-inspired', 'WAOA-inspired']:
    z = np.load(f'{O}logs/diagnostic/base/{p}/777.npz'); c = z['convergence']; ax.plot(range(c.shape[1]), (c/c[:, :1]).mean(axis=0), marker='o', ms=3, label=p)
ax.set_xlabel('Optimizer iteration'); ax.set_ylabel('Best fitness / best of initial population'); ax.legend(fontsize=7, ncol=2); ax.grid(alpha=.25); fig.tight_layout(); fig.savefig(O+'figures/fig2_convergence_seed777.png', dpi=300)
c1 = pd.read_csv(O+'experiment_ctrl.csv'); c1 = c1.drop(columns='protocol').rename(columns={'variant': 'protocol'})
allc = pd.concat([core[core.protocol.isin(['MH-HEED', 'MH-HEED-GSO', 'MH-HEED-PSO'])][['protocol', 'HND']], c1[['protocol', 'HND']]])
names = ['MH-HEED', 'MH-HEED-GREEDY', 'MH-HEED-RS', 'MH-HEED-HC', 'MH-HEED-EXH', 'MH-HEED-GSO', 'MH-HEED-PSO']
mu = [allc[allc.protocol == n].HND.mean() for n in names]; hw = [stats.t.ppf(.975, 19)*allc[allc.protocol == n].HND.std(ddof=1)/np.sqrt(20) for n in names]
fig, ax = plt.subplots(figsize=(7, 3.2)); ax.errorbar(range(7), mu, yerr=hw, fmt='o', capsize=3); ax.set_xticks(range(7)); ax.set_xticklabels(names, rotation=30, ha='right'); ax.set_ylabel('Mean HND (rounds)'); ax.grid(alpha=.25, axis='y'); fig.tight_layout(); fig.savefig(O+'figures/fig3_controls_hnd.png', dpi=300)
print('figures written')
