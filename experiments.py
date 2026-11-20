"""Additional controls requested in the review (see README). Reuses simulator.py unchanged.
New CH-selection controls for the MH-HEED family (same HEED-like candidate pool, same routing, same radio/energy model):
  MH-HEED-GREEDY : deterministic top-K of the pool by the HEED-like score (no objective, no randomness)
  MH-HEED-HC     : hill-climbing (best single swap) on the seven-term objective, started from GREEDY, budget 88 evaluations
  MH-HEED-RS     : random search, best of 88 random K-subsets of the pool under the same objective
  MH-HEED-EXH    : exhaustive evaluation of all K-subsets of the pool (exact minimiser of the objective)
Replicates: task seed = base_seed + 100000*rep; topology/energy streams depend on base_seed only,
selection/optimizer streams additionally depend on rep (rep 0 = original streams)."""
import sys, json, csv, time, itertools, argparse, concurrent.futures
import numpy as np
import simulator as s
np.savez_compressed = lambda *a, **k: None      # experiments do not write per-run logs
CFG = s.CFG
NEW = {'MH-HEED-GREEDY': 21, 'MH-HEED-HC': 22, 'MH-HEED-RS': 23, 'MH-HEED-EXH': 24}
CFG['protocol_stream_ids'].update(NEW)
TERMS = ['D_intra', 'D_route', 'P_meanE', 'P_minE', 'CV_load', 'P_repeat', 'P_relay']
BUDGET = CFG['agents'] * (CFG['iterations'] + 1)      # 88 objective evaluations per election

_orig_streams, _orig_choose = s.streams, s.choose
def streams(seed, protocol):
    base, rep = seed % 100000, seed // 100000
    if rep == 0: return _orig_streams(base, protocol)
    p = CFG['protocol_stream_ids'][protocol]
    return [np.random.default_rng(np.random.SeedSequence(x)) for x in ([base, 0], [base, 1], [base, p*10+2, rep], [base, p*10+3, rep])]

def choose(alive, e, initial, dist, ds, duty, last, r, protocol, c, rr, oo):
    if protocol not in NEW: return _orig_choose(alive, e, initial, dist, ds, duty, last, r, protocol, c, rr, oo)
    ids = np.flatnonzero(alive); k = min(len(ids), max(1, int(np.floor(c['ch_fraction']*len(ids)+.5))))
    local = dist[np.ix_(ids, ids)].sum(axis=1)/max(len(ids)-1, 1)
    score = (e[ids]/initial[ids])/(1+local/max(np.linalg.norm(c['area_m']), 1e-12))
    pool = ids[np.argsort(-score, kind='stable')[:min(len(ids), 2*k)]]
    f = lambda sets: s.objective(np.asarray(sets), alive, dist, ds, e, initial, duty, c, True)
    if protocol == 'MH-HEED-GREEDY': return np.sort(pool[:k]), []
    if protocol == 'MH-HEED-EXH':
        sets = np.array(list(itertools.combinations(pool, k))); return np.sort(sets[f(sets).argmin()]), []
    if protocol == 'MH-HEED-RS':
        sets = pool[np.argsort(oo.random((BUDGET, len(pool))), axis=1)[:, :k]]; sets = np.sort(sets, axis=1)
        return sets[f(sets).argmin()], []
    if protocol == 'MH-HEED-HC':
        cur = np.sort(pool[:k]); fc = float(f(cur[None])[0]); used = 1; rest = np.setdiff1d(pool, cur)
        while used < BUDGET and len(rest) and k:
            cand = np.array([np.sort(np.where(cur == a, b, cur)) for a in cur for b in rest])[:BUDGET-used]
            fv = f(cand); used += len(cand); j = fv.argmin()
            if fv[j] < fc - 1e-15: cur, fc = cand[j], float(fv[j]); rest = np.setdiff1d(pool, cur)
            else: break
        return cur, []
s.streams, s.choose = streams, choose

def run(t):
    meta, task = t; r = s.simulate(task); r.pop('log', None); r.update(meta)
    r['base_seed'], r['rep'] = task[2] % 100000, task[2] // 100000; return r

def build(kind):
    core = CFG['core_seeds']; T = []
    if kind == 'ctrl':
        T = [({'variant': p}, ('ctrl', 'base', sd, p)) for p in NEW for sd in core]
    elif kind == 'ablation':
        w = np.array(CFG['fitness_weights'])
        for i, name in enumerate(TERMS):
            v = w.copy(); v[i] = 0; CFG['stress']['abl_no_'+name] = {'fitness_weights': list(v/v.sum())}
        CFG['stress']['abl_equal'] = {'fitness_weights': [1/7]*7}
        for name in ['abl_no_'+n for n in TERMS] + ['abl_equal']:
            T += [({'variant': name}, ('ablation', name, sd, 'MH-HEED-PSO')) for sd in core]
    elif kind == 'repeat':
        for p in ['MH-HEED', 'MH-HEED-PSO', 'MH-HEED-GSO', 'MH-HEED-RS']:
            T += [({'variant': p}, ('repeat', 'base', sd+100000*rep, p)) for rep in range(1, 10) for sd in core]
    elif kind == 'sens':      # one-at-a-time sensitivity, stress seeds 2001-2010
        S = {'sens_default': {}, 'sens_sink_50_50': {'sink_m': [50, 50]}, 'sens_sink_50_100': {'sink_m': [50, 100]}, 'sens_sink_50_250': {'sink_m': [50, 250]},
             'sens_payload_2000': {'payload_bits': 2000}, 'sens_payload_6000': {'payload_bits': 6000}, 'sens_payload_8000': {'payload_bits': 8000},
             'sens_recluster_10': {'recluster_rounds': 10}, 'sens_recluster_40': {'recluster_rounds': 40}}
        CFG['stress'].update(S)
        for sc in S: T += [({'variant': sc}, ('sens', sc, sd, p)) for sd in CFG['stress_seeds'] for p in ['LEACH', 'MH-HEED', 'MH-HEED-RS', 'MH-HEED-PSO']]
    return T

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('kind', choices=['ctrl', 'ablation', 'repeat', 'sens']); ap.add_argument('--workers', type=int, default=2)
    a = ap.parse_args(); T = build(a.kind); out = s.ROOT/'outputs'/f'experiment_{a.kind}.csv'; t0 = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as ex, out.open('w', newline='') as f:
        w = None
        for i, row in enumerate(ex.map(run, T, chunksize=4), 1):
            if w is None: w = csv.DictWriter(f, fieldnames=list(row)); w.writeheader()
            w.writerow(row); f.flush()
            if i % 40 == 0 or i == len(T): print(f'{a.kind} {i}/{len(T)} {time.time()-t0:.0f}s', flush=True)
