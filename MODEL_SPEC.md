# Model specification — reconstruction v1.0.0
All undefined details were explicitly selected for this reconstruction, rather than attributed to the historical simulator. Consult `simulator.py` for executable semantics.

## Radio and energy accounting
Radio model: TX(l,d)=l(Eelec+eps_fs d²) if d<sqrt(eps_fs/eps_mp), else l(Eelec+eps_mp d⁴). RX(l)=l Eelec; DA(l)=l Eda. Eelec=50 nJ/bit; eps_fs=10 pJ/bit/m²; eps_mp=0.0013 pJ/bit/m⁴; Eda=5 nJ/bit. Each alive node generates one payload per round; a member transmits to nearest CH, a CH aggregates its own payload plus all incoming member and child-CH packets, then transmits one fixed-length packet regardless of subtree size (ideal aggregation). Each child-CH incoming payload is charged RX and DA at its parent. Routes are retained until CH replacement or scheduled reclustering. Every setup: each CH broadcasts 200 bits at max affiliated member distance (0 if none), members listen to all CHs and send 200-bit joins to own CH, each CH sends one 200-bit routing setup packet to parent (or BS), relay CH receives that packet. Sink energy excluded. Centralized global state collection and dissemination costs beyond these declared messages are not modeled; no range cap, collision, interference, or retransmission model.

Round-end accounting: start-alive nodes incur all planned actions, including actions exceeding remaining energy; deduct all planned cost, clip at zero. There is no actual delivery/partial-packet success simulation. This convention is not suitable for hardware packet delivery claims. Death recorded at the end of the first exhausting round; numerical zero tolerance 1e-13 J. State transition: alive set -> elect CH if required -> attach -> route -> setup debit -> data debits -> clip -> deaths. CH death triggers election next round, member death changes receiver load; if no alive node, stop. FND/HND/LND computed directly from individual deaths. If time limit reached with survivors, their death rounds are censored. Model assumes sink has unlimited energy.

Hand-check (no setup, 4000 bits): one member M at 20 m from CH A; A forwards to CH B at 40 m; B forwards to sink at 60 m (all free-space branch). M cost 0.216 mJ = TX(4000,20). A cost 0.504 mJ = RX(4000) 0.200 + DA(8000) 0.040 + TX(4000,40) 0.264. B cost 0.584 mJ = RX(4000) 0.200 + DA(8000) 0.040 + TX(4000,60) 0.344. Total = 1.304 mJ. This is an illustrative route with specified distances, not an observed run or hardware trace.

## CH and routing rules
K=min(alive_count,max(1,floor(0.10*alive_count+0.5))). LEACH-like pool excludes nodes selected in the current 200-round epoch (10 times 20-round interval); if fewer than K remain, use all alive nodes. MH-HEED-like pool ranks alive nodes by (E_i/E_i0)/(1+mean distance to other alive nodes / field diagonal); top min(N_alive,2K) retained, ties stable by node index. Baseline selection is uniform sampling without replacement from the respective pool. This is a centralized simplification rather than the canonical probabilistic LEACH threshold or distributed iterative HEED election. APSO uses LEACH-like pool; WAOA uses HEED-like pool. Members attach nearest selected CH (ties to lower node ID). LEACH family sends CH packets directly to BS. MH family and WAOA use a next hop from other alive CHs strictly closer to sink, which ensures cycle prevention. Candidate cost TX(l,d_ij)*(1+E_initial_nominal/max(E_j,1e-12)); compare to direct TX(l,d_iBS), choose minimum relay if strictly lower, else direct. If node dies at round end, routes with a dead CH cannot continue into next round.

## Objective
Decode top-K unique priorities over selected candidate pool. Let D=max(field diagonal, maximum node-sink distance, 1e-12), C selected CH set, M alive non-CH members; E_i0 initial energy; L_i number of non-CH members assigned to CH i; q_i prior CH appointments; R_i number of immediate child CHs forwarding to CH i.
1. D_intra = sum member distance to nearest CH / max(|M|,1) / D (0 if no members).
2. D_route = mean selected CH next-hop distance / D.
3. P_meanE = 1 - mean(clip(E_i/E_i0,0,1)).
4. P_minE = 1 - min(clip(E_i/E_i0,0,1)).
5. CV_load = population SD(L_i) / max(mean(L_i),1).
6. P_repeat = mean(q_i for selected CHs) / max(max q_all,1).
7. P_relay = max(R_i) / max(K-1,1).
Fitness = weighted sum 0.26,0.14,0.28,0.12,0.08,0.07,0.05 (minimize). Normalization is explicit, not all terms necessarily bounded by one. There is no extra nominal-lifetime factor. The simulator evaluates the initial population plus ten updates: 88 candidate evaluations per election for each optimizer, without historical cache assumptions. No objective term ablation was performed.

## Pseudocode
```
for protocol, scenario, seed in schedule:
  create shared topology and initial energy using dedicated streams
  while alive and round <= max_rounds:
    if interval elapsed OR retained CH dead:
      construct protocol-specific eligible pool, calculate K
      LEACH or MH-HEED: sample K distinct eligible IDs
      LEACH-PSO or MH-HEED-PSO: minimize fitness with PSO priorities
      LEACH-GSO or MH-HEED-GSO: minimize fitness with GSO priorities
      APSO-inspired: PSO plus decaying Gaussian perturbation
      WAOA-inspired: whale-inspired priority update; same MH route rule
      assign members, construct routes, debit setup
    else: keep CH/route; remove dead members and adjust load
    charge planned data TX/RX/DA and update residual energy
    record newly dead nodes, CH, route and energy state
  compute survival metrics, log and runtime
```
PSO: uniform initial priorities in [0,1], velocity zero; inertia 0.9 to 0.4, coefficients 1.8, speed clipped ±0.25, position clipped [0,1]; pbest/gbest improve only on strict decoded fitness. APSO proxy adds decaying Gaussian perturbation SD 0.03*(1-j/10). GSO: initial luciferin=1, update L=0.4L+0.6/(1+F); neighborhood radius initially sqrt(dimension), step 0.12 toward better-luciferin neighbor sampled by positive difference; radius gain 0.08 targeting five neighbors, clipped [0,sqrt(dim)]. Global best retained for final decoded selection. Whale proxy: a=2*(1-j/10), A=2aU-a, C=2U', anchor best if |A|<1 else randomly chosen current agent, x_new=anchor-A*abs(C*anchor-x); clipped [0,1]. No published ant stage reproduced; WAOA-inspired name retained only to match historical comparator category, not to imply original algorithm. All optimization code is a newly specified local approximation.

## Compressed log format
Each NPZ has arrays `start_round,span` (blocks); `energy_before,regular_debit,setup_debit,parent,role,tx_bits,rx_bits,aggregation_bits` (block × node). From start round r, for offset j in [0,span-1], E_before=max(0,E_start-j*regular-(setup if j>0 else 0)), E_after=max(0,E_before-regular-(setup if j=0 else 0)). Apply zero tolerance. Planned debit and actual clipped consumption distinguish energy insufficiency. Role 0 dead, 1 member, 2 CH, 3 relay CH. Parent -1 BS, -2 dead/no route, otherwise zero-based node ID. TX/RX/aggregation bit fields are data bits only; setup joules logged separately. CH, relay, and energy are recoverable for **every round** without losing state. A block ends on earliest node death or next scheduled election; fixed routing and cost within each block. Arrays `xy,initial,death,ch_duty` hold node-level values; `convergence_round,convergence` record election-round fitness sequences. Export script writes per-node per-round CSV. Runtime includes simulation and writing compressed logs, not a measurement of on-node computation.
